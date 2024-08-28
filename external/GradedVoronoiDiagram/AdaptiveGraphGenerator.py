from drbutil import *
from scipy.spatial import Voronoi
from scipy.spatial import Delaunay
from scipy.spatial import KDTree

try:
    import igl
except ImportError:
    iglFound = False
else:
    iglFound = True

try:
    from embreex import rtcore as rtc
    from embreex import rtcore_scene as rtcs
    from embreex.mesh_construction import TriangleMesh
except ImportError:
    embreeFound = False
else:
    embreeFound = True


class AdaptiveGraphGenerator:

    def __init__(self, filePath):

        self.verts, self.ts, flxIdxs, forceVecs, fixIdxs, self.vmStress, pStress, pStressE, sMats = loadStressFile(filePath)
        self.scalars = icdf(self.vmStress) #normZeroToOne(self.vmStress)
        self.filePath = filePath.replace('.stress', '')

        self.nDim = self.verts.shape[1]

        self.boundaryVertMask = np.zeros(len(self.verts), np.bool_)
        if self.nDim == 2:
            self.edges = facesToEdges(self.ts, False)
            unq, inv, cnt = np.unique(cantorPiKV(self.edges), return_inverse = True, return_counts = True)
            self.boundaryEdgeMask = cnt[inv] == 1
            self.boundaryVertMask[np.unique(self.edges[self.boundaryEdgeMask].ravel())] = True
        else:
            self.tris = tetrasToFaces(self.ts)
            self.edges = facesToEdges(self.tris, True)
            unq, inv, cnt = np.unique(cantorPiKV(self.tris), return_inverse = True, return_counts = True)
            self.boundaryTrisMask = cnt[inv] == 1
            self.boundaryVertMask[np.unique(self.tris[self.boundaryTrisMask].ravel())] = True

        self.mnmx = np.float32([self.verts.min(axis=0), self.verts.max(axis=0)])
        self.BBdiag = norm(self.mnmx[1] - self.mnmx[0])
        self.sf = self.BBdiag / 1000
        self.mel = norm(self.verts[self.edges[:,0]] - self.verts[self.edges[:,1]]).mean()
        
    def createSamples(self, scale = 0.05, ratio = 0.25):
        self.scale = scale
        self.ratio = 1 if np.allclose(self.scalars, 1) else np.clip(ratio, 0.001, 1)
        sMax = self.BBdiag * self.scale
        self.s = self.scalars * (sMax * self.ratio - sMax) + sMax

        # auxiliary samples around the object       
        auxSamples, _ = (generatePointsOnCircle(12), None) if self.nDim == 2 else generateIcoSphere(1, True)
        self.samples = auxSamples * self.BBdiag + self.mnmx.mean(axis=0)
        self.sRadii = np.zeros(len(auxSamples), np.float32)

        # hull vertices from the input
        for sample, r in zip(self.verts[self.boundaryVertMask], self.s[self.boundaryVertMask]):
            dsts = norm(self.samples - sample)
            if not len(dsts) or np.all(dsts > r) and np.all(dsts > self.sRadii):
                self.samples = np.vstack([self.samples, sample])
                self.sRadii = np.concatenate([self.sRadii, [r]])

        trysLeft = 5
        nSamplesPerRound = int(len(self.ts) * self.mel / sMax)
        nMin = nSamplesPerRound/1000
        while trysLeft:
            tsSubset = self.ts[np.random.randint(0, len(self.ts), nSamplesPerRound)]
            
            bWeights = np.random.rand(nSamplesPerRound, self.nDim + 1)**self.nDim
            bWeights /= bWeights.sum(axis=1).reshape(-1,1)

            samples = (bWeights[:,:,np.newaxis] * self.verts[tsSubset]).sum(axis=1)
            radii = (bWeights * self.s[tsSubset]).sum(axis=1)

            kdtAll = KDTree(self.samples)
            kdtNew = KDTree(samples)

            nnssAll = kdtAll.query_ball_point(samples, radii)
            nnssNew = kdtNew.query_ball_point(samples, radii)
            nnssNew = list(map(set, nnssNew))

            sIdxsToAdd = set([])
            for sIdx in tqdm(range(nSamplesPerRound), total = nSamplesPerRound, ascii=True, desc='sampling'):               
                if not len(nnssAll[sIdx]) and nnssNew[sIdx].isdisjoint(sIdxsToAdd):
                    sIdxsToAdd.add(sIdx)

            sIdxs = np.int32(list(sIdxsToAdd))
            self.samples = np.vstack([self.samples, samples[sIdxs]])
            self.sRadii = np.concatenate([self.sRadii, radii[sIdxs]])

            if len(sIdxsToAdd) < nSamplesPerRound/1000:
                trysLeft -= 1

            print('nSamples:', len(self.samples), 'nAdded (min):', len(sIdxsToAdd), '(%d)'%nMin, 'trysLeft:', trysLeft)

    def computeVoronoi(self):
        self.vr = Voronoi(self.samples)
        self.vrVerts = self.vr.vertices       
        
        cells = [self.vr.regions[rIdx] for rIdx in self.vr.point_region] if self.nDim == 2 else self.vr.ridge_vertices
        self.cells = [np.int64(cell) for cell in cells if not -1 in cell]
        cellEdges = [faceToEdges(cell) for cell in self.cells]

        vs = self.verts[self.ts]
        if self.nDim == 2:
            boundaryEdgeVerts = self.verts[self.edges[self.boundaryEdgeMask]]
            vrIo = np.bool_([pointInTriangles2D(vs, p) for p in self.vrVerts])
        else:
            boundaryTris = self.tris[self.boundaryTrisMask]
            boundaryTrisVerts = self.verts[boundaryTris]
            
            if iglFound:
                vrIo = np.abs(igl.fast_winding_number_for_meshes(np.float64(self.verts), boundaryTris, self.vrVerts)) > 0.5
            else:
                subTVerts = vs[:,[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]].reshape(-1, self.nDim, self.nDim)
                vrIo = np.bool_([np.all(simpleDets3x3(subTVerts - vert.reshape(1,3)).reshape(-1,4) >= 0, axis=1).any() for vert in self.vrVerts])

            edges = filterForUniqueEdges(np.vstack(cellEdges))
            edges = edges[vrIo[edges].sum(axis=1) == 1]
            origins = self.vrVerts[edges[:,0]]
            dirs = normVec(self.vrVerts[edges[:,1]] - origins)

            if embreeFound:
                embreeDevice = rtc.EmbreeDevice()
                scene = rtcs.EmbreeScene(embreeDevice)
                mesh = TriangleMesh(scene, boundaryTrisVerts)
                dsts = scene.run(np.float32(origins - dirs * eps), np.float32(dirs), query='DISTANCE')
            else:
                dsts = [min(intersectTrianglesWithRay(boundaryTrisVerts, o - d * eps, d)) for e, o, d in zip(edges, origins, dirs)]

            eDsts = {cantorPi(e[0], e[1]): dst for e, dst in zip(edges, dsts)}

        vrEdges = []
        cutEdgeVertIdxs = {}
        newVerts = []
        
        for cIdx, cell in enumerate(self.cells):
            newEdge = []
            for e in cellEdges[cIdx]:
                pIo, qIo = vrIo[e]
                    
                if pIo and qIo:
                    vrEdges.append(e)
                elif pIo or qIo:
                    eHash = cantorPi(e[0], e[1])

                    p, q = self.vrVerts[e]
                    d = normVec(q-p)
                    newVert = intersectEdgesWithRay2D(boundaryEdgeVerts, p, d) if self.nDim == 2 else (p + d * eDsts[eHash])

                    if eHash in cutEdgeVertIdxs.keys():
                        newVertIdx = cutEdgeVertIdxs[eHash]
                    else:
                        newVertIdx = len(self.vrVerts) + len(newVerts)
                        cutEdgeVertIdxs[eHash] = newVertIdx
                    
                    newVerts.append(newVert)
                    newEdge.append(newVertIdx)

                    vrEdges.append([e[[pIo,qIo]][0], newVertIdx])

            if len(newEdge) == 2:
                vrEdges.append(newEdge)

        self.vrVerts = np.vstack([self.vrVerts, newVerts])[np.unique(vrEdges)]
        self.vrEdges = reIndexIndices(np.int32(vrEdges))

    def computeDelaunay(self):
        dl = Delaunay(self.samples)

        vs = self.verts[self.ts]
        centroids = dl.points[dl.simplices].mean(axis=1)
        if self.nDim == 2:
            dlIo = np.bool_([pointInTriangles2D(vs, p) for p in centroids])
        else:
            if iglFound:
                self.ctrs = centroids
                self.wns = igl.fast_winding_number_for_meshes(np.float64(self.verts), self.tris[self.boundaryTrisMask], centroids)
                dlIo = np.abs(igl.fast_winding_number_for_meshes(np.float64(self.verts), self.tris[self.boundaryTrisMask], centroids)) > 0.5
            else:
                subTVerts = vs[:,[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]].reshape(-1, self.nDim, self.nDim)
                dlIo = np.bool_([np.all(simpleDets3x3(subTVerts - vert.reshape(1,3)).reshape(-1,4) >= 0, axis=1).any() for vert in centroids])
        
        dlEdges = facesToEdges(dl.simplices[dlIo]) if self.nDim == 2 else tetsToEdges(dl.simplices[dlIo])
        
        self.dlVerts = dl.points[np.unique(dlEdges)]
        self.dlEdges = reIndexIndices(dlEdges)
        self.dlCells = reIndexIndices(dl.simplices[dlIo])

    def exportResults(self, filePath = None):
        filePath = self.filePath if filePath is None else filePath
        dlVerts = pad2Dto3D(self.dlVerts) if self.nDim == 2 else self.dlVerts
        writeObjFile(filePath + "_Delaunay.obj", dlVerts, [], self.dlEdges)
        vrVerts = pad2Dto3D(self.vrVerts) if self.nDim == 2 else self.vrVerts
        writeObjFile(filePath + "_Voronoi.obj", vrVerts, [], self.vrEdges)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', type=str)
    parser.add_argument('s', type=float)
    parser.add_argument('r', type=float)
    args = parser.parse_args()
    if not os.path.exists(args.inputFile) or not args.inputFile.endswith('.stress'):
        parser.error("Input file does not exist or not in correct format.")

    agg = AdaptiveGraphGenerator(args.inputFile)
    st = time()
    agg.createSamples(args.s, args.r)
    print('t Sampling:\t', time() - st)
    st = time()
    agg.computeVoronoi()
    print('t Voronoi:\t', time() - st)
    st = time()
    agg.computeDelaunay()
    print('t Delaunay:\t', time() - st)
    agg.exportResults()
