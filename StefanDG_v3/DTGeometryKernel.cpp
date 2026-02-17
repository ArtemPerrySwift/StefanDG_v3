#include "DTGeometryKernel.h"
#include "GMSHProxy.h"

#include <unordered_map>
#include <stdexcept>
#include <numeric>

extern "C" {
#include "gmshc.h"
}

namespace DTGeometryKernel
{
    inline static void throwLastError()
    {
        int ierr = 0;
        char* api_error_;
        gmshLoggerGetLastError(&api_error_, &ierr);
        if (ierr) throw std::runtime_error("Could not get last error");
        std::string error = std::string(api_error_);
        gmshFree(api_error_);
        throw std::runtime_error(error);
    }

        void assignVolumesMaterialsPhases(const int volumesTags[], const size_t nVolumes, const MaterialPhase* materialPhaseIt, const uint8_t nMaterialsPhases, const MaterialPhase* volumesMaterialsPhasesPtrs[])
        {
            for (uint8_t i = 0; i < nMaterialsPhases; ++i)
            {
                int* dimTags;
                size_t nDimTags;

                int ierr = 0;
                gmshModelGetEntitiesForPhysicalName(materialPhaseIt->name, &dimTags, &nDimTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                size_t nMaterialPhaseEntities = nDimTags / 2;
                for (size_t j = 0; j < nMaterialPhaseEntities; ++j)
                {
                    ++dimTags;
                    int tag = *dimTags;
                    const int* volumeTagsIt = volumesTags;

                    size_t k = 0;
                    while (tag != *volumeTagsIt)
                    {
                        ++k;
                        ++volumeTagsIt;
                    }

                    volumesMaterialsPhasesPtrs[k] = materialPhaseIt;
                    ++dimTags;
                }

                ++materialPhaseIt;
            }
        }

        size_t countNonconformalVolumesInterfaces()
        {
            int* dimTags;
            size_t nDimTags;
            GMSHProxy::model::getPhysicalGroups(dimTags, nDimTags, constants::DIMENSION_2D);

            size_t nNonconformalInterfaces = 0;
            int* dimTagIt = dimTags;
            size_t nSurfaces = nDimTags / 2;
            for (size_t i = 0; i < nSurfaces; ++i)
            {
                ++dimTagIt;
                char* physicalName;
                GMSHProxy::model::getPhysicalName(constants::DIMENSION_2D, *dimTagIt, physicalName);
                if (*physicalName == 'C')
                {
                    ++nNonconformalInterfaces;
                }

                GMSHProxy::free(physicalName);
                ++dimTagIt;
            }
            GMSHProxy::free(dimTags);

            return nNonconformalInterfaces;
        }

        void extractFaceNodesTags(const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES], const uint8_t faceLocalIndex, size_t* faceNodesTagsIt)
        {
            switch (faceLocalIndex)
            {
            case 0:
            {
                *faceNodesTagsIt = tetrahedronNodesTags[0];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[2];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[1];
                ++faceNodesTagsIt;

                break;
            }
            case 1:
            {
                *faceNodesTagsIt = tetrahedronNodesTags[0];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[1];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[3];
                ++faceNodesTagsIt;

                break;
            }
            case 2:
            {
                *faceNodesTagsIt = tetrahedronNodesTags[0];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[3];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[2];
                ++faceNodesTagsIt;

                break;
            }
            case 3:
            {
                *faceNodesTagsIt = tetrahedronNodesTags[3];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[1];
                ++faceNodesTagsIt;

                *faceNodesTagsIt = tetrahedronNodesTags[2];
                ++faceNodesTagsIt;

                break;
            }
            default:
                throw std::runtime_error("Wrong face index");
            }
        }

        void extractBoundaryTetrahedrons(const double tetrahedronsNativeJacobians[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         const Coordinates tetrahedronsInitPoints[],
                                         const size_t* boundaryTetrahedronsIndexIt,
                                         const size_t nBoundaryTetrahedrons,
                                         double (*boundaryTetrahedronsNativeJacobiansIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         Coordinates* boundaryTetrahedronInitPointsIt)
        {
            for (size_t i = 0; i < nBoundaryTetrahedrons; ++i)
            {
                std::copy_n(tetrahedronsNativeJacobians[*boundaryTetrahedronsIndexIt],
                    LocalCoordinates3D::COUNT * Coordinates::COUNT,
                    *boundaryTetrahedronsNativeJacobiansIt);

                *boundaryTetrahedronInitPointsIt = tetrahedronsInitPoints[*boundaryTetrahedronsIndexIt];

                ++boundaryTetrahedronInitPointsIt;
                ++boundaryTetrahedronsNativeJacobiansIt;
                ++boundaryTetrahedronsIndexIt;
            }
        }

        static void freeArrays(int** arrOfArrIt, size_t nArr)
        {
            for (size_t i = 0; i < nArr; ++i)
            {
                gmshFree(*arrOfArrIt);
                ++arrOfArrIt;
            }
        }

        static void translate(const size_t* templateMeshTagIt, const int* templateModelTagIt, const size_t nTemplateTags, const size_t* meshTagIt, int* modelTagIt, const size_t nOutTags)
        {
            std::unordered_map<size_t, int> modelTagByMeshTag;
            modelTagByMeshTag.reserve(nTemplateTags);
            for (size_t i = 0; i < nTemplateTags; ++i)
            {
                modelTagByMeshTag.emplace(*templateMeshTagIt, *templateModelTagIt);

                ++templateMeshTagIt;
                ++templateModelTagIt;
            }

            for (size_t i = 0; i < nOutTags; ++i)
            {
                *modelTagIt = modelTagByMeshTag.at(*meshTagIt);

                ++modelTagIt;
                ++meshTagIt;
            }
        }

        void addTriangleEntities(const size_t* trianglesEdgeTagIt, const size_t nTriangles, const int* triagnlesEdgesNodeTagIt, int* outEntitiesDimTagIt)
        {
            int ierr;
            int elementEdgesTags[constants::triangle::N_EDGES];
            std::unordered_map<size_t, int> lineTagByEdgeTag;
            lineTagByEdgeTag.reserve(nTriangles);

            for (size_t i = 0; i < nTriangles; ++i)
            {
                for (uint8_t j = 0; j < constants::triangle::N_EDGES; ++j)
                {
                    size_t edgeTag = *trianglesEdgeTagIt;
                    auto it = lineTagByEdgeTag.find(edgeTag);
                    if (it == lineTagByEdgeTag.end())
                    {
                        int lineTag = gmshModelOccAddLine(*triagnlesEdgesNodeTagIt, *(triagnlesEdgesNodeTagIt + 1), GMSHProxy::model::DEFAULT_TAG, &ierr);
                        if (ierr)
                        {
                            throwLastError();
                        }

                        lineTagByEdgeTag.emplace(edgeTag, lineTag);
                        elementEdgesTags[j] = lineTag;
                    }
                    else
                    {
                        elementEdgesTags[j] = it->second;
                    }

                    ++trianglesEdgeTagIt;
                    triagnlesEdgesNodeTagIt += 2;
                }

                int curveLoopTag = gmshModelOccAddCurveLoop(elementEdgesTags, constants::triangle::N_EDGES, GMSHProxy::model::DEFAULT_TAG, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                *outEntitiesDimTagIt = constants::DIMENSION_2D;
                ++outEntitiesDimTagIt;
                *outEntitiesDimTagIt = gmshModelOccAddPlaneSurface(&curveLoopTag, 1, GMSHProxy::model::DEFAULT_TAG, &ierr);
                ++outEntitiesDimTagIt;

                if (ierr)
                {
                    throwLastError();
                }
            }
        }

        static void buildEntitiesBasedOnSurfaceTriangles(const int surfaceTag, const char sourceModelName[], const char reciverModelName[], int*& outEntitiesDimTags, size_t& nOutEntitiesDimTag)
        {
            int ierr;
            size_t* nodesMeshTags;
            double* nodesCoordinates, * parametricCoordinates;
            size_t nNodes, nNodesCoordinates, nParametricCoordinates;
            gmshModelMeshGetNodes(&nodesMeshTags,
                &nNodes,
                &nodesCoordinates,
                &nNodesCoordinates,
                &parametricCoordinates,
                &nParametricCoordinates,
                constants::DIMENSION_2D,
                surfaceTag,
                GMSHProxy::TRUE_C,
                GMSHProxy::FALSE_C,
                &ierr);


            if (ierr)
            {
                throwLastError();
            }

            gmshModelSetCurrent(reciverModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            int* nodesModelTags = new int[nNodes];
            GMSHProxy::model::addNodes(nodesCoordinates, nNodes, nodesModelTags);

            gmshModelSetCurrent(sourceModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshFree(nodesCoordinates);
            gmshFree(parametricCoordinates);

            int surfaceDimTag[] = { constants::DIMENSION_2D,  surfaceTag };
            gmshModelMeshCreateEdges(surfaceDimTag, 2, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            size_t* edgesNodesMeshTags;
            size_t nEdgesNodesTags;
            gmshModelMeshGetElementEdgeNodes(GMSHProxy::model::mesh::TRIANGLE_TYPE, &edgesNodesMeshTags, &nEdgesNodesTags, surfaceTag, GMSHProxy::FALSE_C, GMSHProxy::TASK_DEFAULT, GMSHProxy::N_TASK_DEFAULT, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            size_t* edgesTags;
            int* edgesOrientations;
            size_t nEdges, nOrientations;

            gmshModelMeshGetEdges(edgesNodesMeshTags, nEdgesNodesTags, &edgesTags, &nEdges, &edgesOrientations, &nOrientations, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshFree(edgesOrientations);

            int* edgesNodesModelTags = new int[nEdgesNodesTags];
            translate(nodesMeshTags, nodesModelTags, nNodes, edgesNodesMeshTags, edgesNodesModelTags, nEdgesNodesTags);

            gmshFree(nodesMeshTags);
            delete[] nodesModelTags;

            const size_t nTriangles = nEdges / constants::triangle::N_EDGES;
            nOutEntitiesDimTag = 2 * nTriangles;

            gmshModelSetCurrent(reciverModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            outEntitiesDimTags = new int[nOutEntitiesDimTag];

            addTriangleEntities(edgesTags, nTriangles, edgesNodesModelTags, outEntitiesDimTags);

            delete[] edgesNodesModelTags;
            gmshFree(edgesTags);

            gmshModelSetCurrent(sourceModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }
        }

        void determineFragmentsTrianglesIndexes(const int* fragmentsTagsIt,
            const size_t nFragments,
            const int* const* trianglesFragmentsDimTagIt,
            size_t* nTrianglesFragmentsDimTagsIt,
            const size_t nObjectTriangles,
            const size_t nToolTriangles,
            size_t(*&fragmentsTrianglesIndexes)[N_CROSS_ELEMENTS])
        {
            std::unordered_map<int, size_t> fragmentIndexByFragmentTag;
            fragmentIndexByFragmentTag.reserve(nFragments);

            for (size_t i = 0; i < nFragments; ++i)
            {
                fragmentIndexByFragmentTag[*fragmentsTagsIt] = i;
                ++fragmentsTagsIt;
            }

            for (size_t i = 0; i < nObjectTriangles; ++i)
            {
                size_t nTriangleFragments = *nTrianglesFragmentsDimTagsIt / 2;
                const int* triangleFragmentsDimTagIt = *trianglesFragmentsDimTagIt;

                for (size_t j = 0; j < nTriangleFragments; ++j)
                {
                    ++triangleFragmentsDimTagIt;
                    size_t fragmentIndex = fragmentIndexByFragmentTag.at(*triangleFragmentsDimTagIt);
                    fragmentsTrianglesIndexes[fragmentIndex][0] = i;

                    ++triangleFragmentsDimTagIt;
                }

                ++trianglesFragmentsDimTagIt;
                ++nTrianglesFragmentsDimTagsIt;
            }


            for (size_t i = 0; i < nToolTriangles; ++i)
            {
                size_t nTriangleFragments = *nTrianglesFragmentsDimTagsIt / 2;
                const int* triangleFragmentsDimTagIt = *trianglesFragmentsDimTagIt;

                for (size_t j = 0; j < nTriangleFragments; ++j)
                {
                    ++triangleFragmentsDimTagIt;
                    size_t fragmentIndex = fragmentIndexByFragmentTag.at(*triangleFragmentsDimTagIt);
                    fragmentsTrianglesIndexes[fragmentIndex][1] = i;

                    ++triangleFragmentsDimTagIt;
                }

                ++trianglesFragmentsDimTagIt;
                ++nTrianglesFragmentsDimTagsIt;
            }
        }

        void buildFragmentsModel(const int surafcesTags[N_CROSS_ELEMENTS],
            const char newModelName[],
            int*& fragmentTags,
            size_t& nFragments,
            size_t(*&surfacesFragmentsTrianglesIndexes)[N_CROSS_ELEMENTS])
        {
            char* mainModelName;
            int ierr;
            gmshModelGetCurrent(&mainModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshModelAdd(newModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshModelSetCurrent(mainModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            int* objectTrianglesEntitiesDimTags;
            size_t nObjectTrianglesEntitiesDimTags;
            int* toolTrianglesEntitiesDimTags;
            size_t nToolTrianglesEntitiesDimTags;

            buildEntitiesBasedOnSurfaceTriangles(surafcesTags[0], mainModelName, newModelName, objectTrianglesEntitiesDimTags, nObjectTrianglesEntitiesDimTags);
            buildEntitiesBasedOnSurfaceTriangles(surafcesTags[1], mainModelName, newModelName, toolTrianglesEntitiesDimTags, nToolTrianglesEntitiesDimTags);

            int* fragmentsDimTags;
            size_t nFragemntsDimTags;
            int** fragmentsDimTagsMap;
            size_t* nFragmentsDimTagsMap;
            size_t nnFragemntsDimTagsMap;

            gmshModelSetCurrent(newModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshModelOccFragment(objectTrianglesEntitiesDimTags, nObjectTrianglesEntitiesDimTags,
                toolTrianglesEntitiesDimTags, nToolTrianglesEntitiesDimTags,
                &fragmentsDimTags, &nFragemntsDimTags,
                &fragmentsDimTagsMap, &nFragmentsDimTagsMap, &nnFragemntsDimTagsMap,
                GMSHProxy::model::DEFAULT_TAG,
                GMSHProxy::TRUE_C,
                GMSHProxy::TRUE_C,
                &ierr);

            if (ierr)
            {
                throwLastError();
            }

            delete[] objectTrianglesEntitiesDimTags;
            delete[] toolTrianglesEntitiesDimTags;

            nFragments = nFragemntsDimTags / 2;
            fragmentTags = (int*)gmshMalloc(sizeof(int) * nFragments);
            GMSHProxy::model::extractTags(fragmentsDimTags, nFragments, fragmentTags);
            gmshFree(fragmentsDimTags);

            surfacesFragmentsTrianglesIndexes = (size_t(*)[N_CROSS_ELEMENTS])gmshMalloc(sizeof(size_t) * nFragments * N_CROSS_ELEMENTS);

            determineFragmentsTrianglesIndexes(fragmentTags, nFragments, fragmentsDimTagsMap, nFragmentsDimTagsMap, nObjectTrianglesEntitiesDimTags / 2, nToolTrianglesEntitiesDimTags / 2, surfacesFragmentsTrianglesIndexes);

            freeArrays(fragmentsDimTagsMap, nnFragemntsDimTagsMap);
            gmshFree(fragmentsDimTagsMap);
            gmshFree(nFragmentsDimTagsMap);

            gmshModelOccSynchronize(&ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshModelSetCurrent(mainModelName, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            gmshFree(mainModelName);

        }

        void determineFragmentsTetrahderonsIndexes(const size_t nFragments,
                                                   const size_t(*fragmentsFacesIndexesIt)[N_CROSS_ELEMENTS],
                                                   const size_t* const intefaceTetrahedronsIndexes[N_CROSS_ELEMENTS],
                                                   size_t(*fragmentsTetrahedronsIndexesIt)[N_CROSS_ELEMENTS])
        {
            for (size_t i = 0; i < nFragments; ++i)
            {
                const size_t* fragmentFacesIndexes = *fragmentsFacesIndexesIt;
                size_t* fragmentTetrahedronsIndexes = *fragmentsTetrahedronsIndexesIt;

                fragmentTetrahedronsIndexes[0] = intefaceTetrahedronsIndexes[0][fragmentFacesIndexes[0]];
                fragmentTetrahedronsIndexes[1] = intefaceTetrahedronsIndexes[1][fragmentFacesIndexes[1]];

                ++fragmentsFacesIndexesIt;
                ++fragmentsTetrahedronsIndexesIt;
            }
        }


        static void determineFaceTetrahedronConnections(const size_t* tetrahedronsFacesTags,
            const size_t nTetrahedrons,
            const size_t* const* boundariesFacesTags,
            const size_t* nBoundariesFaces,
            const size_t nBoundaries,
            const size_t nInteriorFaces,
            const size_t nFaces,
            size_t(*interiorFacesTetrahedronsIndexes)[2],
            uint8_t(*interiorFacesLocalIndexes)[2],
            size_t* boundariesFacesTetrahedronsIndexes[],
            uint8_t* boundariesFacesLocalIndexes[])
        {
            std::unordered_map<size_t, size_t> faceIndexeByFaceTag;
            faceIndexeByFaceTag.reserve(nFaces);

            for (size_t iTetrahedron = 0, iFace = 0; iTetrahedron < nTetrahedrons; ++iTetrahedron)
            {
                for (uint8_t iFaceLocal = 0; iFaceLocal < constants::tetrahedron::N_FACES; ++iFaceLocal, ++iFace)
                {
                    auto faceIndexesByFaceTagIt = faceIndexeByFaceTag.find(*tetrahedronsFacesTags);
                    if (faceIndexesByFaceTagIt == faceIndexeByFaceTag.end())
                    {
                        faceIndexeByFaceTag.emplace(*tetrahedronsFacesTags, iFace);
                    }
                    else
                    {
                        size_t* interiorFaceTetrahedronsIndexes = *interiorFacesTetrahedronsIndexes;
                        size_t previousFaceIndex = faceIndexeByFaceTag.at(*tetrahedronsFacesTags);
                        interiorFaceTetrahedronsIndexes[0] = previousFaceIndex / constants::tetrahedron::N_FACES;
                        interiorFaceTetrahedronsIndexes[1] = iTetrahedron;

                        uint8_t* interiorFaceLocalIndexes = *interiorFacesLocalIndexes;
                        interiorFaceLocalIndexes[0] = previousFaceIndex % constants::tetrahedron::N_FACES;
                        interiorFaceLocalIndexes[1] = iFaceLocal;

                        ++interiorFacesTetrahedronsIndexes;
                        ++interiorFacesLocalIndexes;
                    }
                    ++tetrahedronsFacesTags;
                }
            }

            for (size_t iBoundary = 0; iBoundary < nBoundaries; ++iBoundary)
            {
                const size_t* boundaryFacesTagsIt = *boundariesFacesTags;
                size_t* boundaryFacesTetrahedronsIndexesIt = *boundariesFacesTetrahedronsIndexes;
                uint8_t* boundaryFacesLocalIndexesIt = *boundariesFacesLocalIndexes;
                const size_t nBoundaryFaces = *nBoundariesFaces;

                for (size_t iFace = 0; iFace < nBoundaryFaces; ++iFace)
                {
                    size_t faceIndex = faceIndexeByFaceTag.at(*boundaryFacesTagsIt);
                    *boundaryFacesTetrahedronsIndexesIt = faceIndex / constants::tetrahedron::N_FACES;
                    *boundaryFacesLocalIndexesIt = faceIndex % constants::tetrahedron::N_FACES;

                    ++boundaryFacesTagsIt;
                    ++boundaryFacesTetrahedronsIndexesIt;
                    ++boundaryFacesLocalIndexesIt;
                }

                ++boundariesFacesTags;
                ++boundariesFacesTetrahedronsIndexes;
                ++boundariesFacesLocalIndexes;
                ++nBoundariesFaces;
            }
        }

        static void getBoundariesFacesTags(const int* boundariesTagsIt,
            const size_t nBoundaries,
            size_t** boundariesFacesTagsIt,
            size_t* nBoundariesFacesIt)
        {
            size_t* elementTags;
            size_t* nodeTags;
            size_t nElements;

            for (size_t i = 0; i < nBoundaries; ++i)
            {
                GMSHProxy::model::mesh::getElementsByType(GMSHProxy::model::mesh::TRIANGLE_TYPE, elementTags, nElements, nodeTags, *boundariesTagsIt);
                GMSHProxy::free(elementTags);

                GMSHProxy::model::mesh::getFaces(GMSHProxy::model::mesh::TRIANGLE_FACE_TYPE, nodeTags, nElements * constants::triangle::N_NODES, boundariesFacesTagsIt);
                GMSHProxy::free(nodeTags);

                *nBoundariesFacesIt = nElements;

                ++boundariesTagsIt;
                ++boundariesFacesTagsIt;
                ++nBoundariesFacesIt;
            }
        }

        static void deleteBoundariesFacesTags(const size_t nBoundaries, size_t** boundariesFacesTagsIt)
        {
            for (size_t i = 0; i < nBoundaries; ++i)
            {
                GMSHProxy::free(*boundariesFacesTagsIt);
                ++boundariesFacesTagsIt;
            }
        }

        static void allocateMemory(const size_t* nBoundariesFacesIt,
            const size_t nBoundaries,
            size_t** boundariesFacesTetrahedronsIndexesIt,
            uint8_t** boundariesFacesLocalIndexesIt)
        {
            for (size_t i = 0; i < nBoundaries; ++i)
            {
                *boundariesFacesTetrahedronsIndexesIt = new size_t[*nBoundariesFacesIt];
                *boundariesFacesLocalIndexesIt = new uint8_t[*nBoundariesFacesIt];

                ++boundariesFacesTetrahedronsIndexesIt;
                ++boundariesFacesLocalIndexesIt;
                ++nBoundariesFacesIt;
            }
        }

        static int createInteriorFacesEntity(const size_t tetrahedronsFacesNodesTags[][constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
            const size_t(*facesTetrahedronsIndexesIt)[2],
            const uint8_t(*faceLocalIndexesIt)[2],
            const size_t nFaces)
        {
            int facesEntityTag = GMSHProxy::model::addDescreteEntity(constants::DIMENSION_2D);
            size_t nFacesNodes = nFaces * constants::triangle::N_NODES;
            size_t* facesNodes = new size_t[nFacesNodes];

            size_t* faceNodesIt = facesNodes;
            for (size_t i = 0; i < nFaces; ++i)
            {
                faceNodesIt = std::copy_n(tetrahedronsFacesNodesTags[**facesTetrahedronsIndexesIt][**faceLocalIndexesIt], constants::triangle::N_NODES, faceNodesIt);
                ++facesTetrahedronsIndexesIt;
                ++faceLocalIndexesIt;
            }

            GMSHProxy::model::mesh::addElementsByType(facesEntityTag, GMSHProxy::model::mesh::TRIANGLE_TYPE, facesNodes, nFacesNodes);
            delete[] facesNodes;

            return facesEntityTag;
        }

        static int createFacesEntity(const size_t tetrahedronsFacesNodesTags[][constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
            const size_t* facesTetrahedronsIndexesIt,
            const uint8_t* faceLocalIndexesIt,
            const size_t nFaces)
        {
            int facesEntityTag = GMSHProxy::model::addDescreteEntity(constants::DIMENSION_2D);
            size_t nFacesNodes = nFaces * constants::triangle::N_NODES;
            size_t* facesNodes = new size_t[nFacesNodes];

            size_t* faceNodesIt = facesNodes;
            for (size_t i = 0; i < nFaces; ++i)
            {
                faceNodesIt = std::copy_n(tetrahedronsFacesNodesTags[*facesTetrahedronsIndexesIt][*faceLocalIndexesIt], constants::triangle::N_NODES, faceNodesIt);
                ++facesTetrahedronsIndexesIt;
                ++faceLocalIndexesIt;
            }

            GMSHProxy::model::mesh::addElementsByType(facesEntityTag, GMSHProxy::model::mesh::TRIANGLE_TYPE, facesNodes, nFacesNodes);
            delete[] facesNodes;

            return facesEntityTag;
        }

        static void createBoundariesFacesEntities(const size_t tetrahedronsFacesNodesTags[][constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
            const size_t* const* boundariesFacesTetrahedronsIndexesIt,
            const uint8_t* const* boundariesFacesLocalIndexesIt,
            const size_t* nBoundariesFacesIt,
            const size_t nBoundaries,
            int* boundariesFacesEntitiesTagsIt)
        {
            for (size_t iBoundary = 0; iBoundary < nBoundaries; ++iBoundary)
            {
                *boundariesFacesEntitiesTagsIt = createFacesEntity(tetrahedronsFacesNodesTags, *boundariesFacesTetrahedronsIndexesIt, *boundariesFacesLocalIndexesIt, *nBoundariesFacesIt);

                ++boundariesFacesTetrahedronsIndexesIt;
                ++nBoundariesFacesIt;
                ++boundariesFacesEntitiesTagsIt;
                ++boundariesFacesLocalIndexesIt;
            }
        }



        static void createBoundariesFacesEntities(const size_t tetrahedronsFacesNodesTags[][constants::tetrahedron::N_FACES][constants::triangle::N_NODES],
                                                  const size_t* const* boundariesFacesTetrahedronsIndexesIt,
                                                  const uint8_t* const* boundariesFacesLocalIndexesIt,
                                                  const size_t* nBoundariesFacesIt,
                                                  const size_t nBoundaries,
                                                  const BoundaryCondition* boundaryConditionIt,
                                                  int* boundariesFacesEntitiesTagsIt)
        {
            for (size_t i = 0; i < nBoundaries; ++i)
            {
                if (*boundaryConditionIt != BoundaryCondition::HOMOGENEOUS_NEWMAN || *boundaryConditionIt != BoundaryCondition::NONCONFORM_INTERAFACE)
                {
                    *boundariesFacesEntitiesTagsIt = createFacesEntity(tetrahedronsFacesNodesTags, *boundariesFacesTetrahedronsIndexesIt, *boundariesFacesLocalIndexesIt, *nBoundariesFacesIt);
                }
                else
                {
                    *boundariesFacesEntitiesTagsIt = GMSHProxy::model::UNSUPPORTED_TAG;
                }


                ++boundariesFacesTetrahedronsIndexesIt;
                ++boundariesFacesLocalIndexesIt;
                ++nBoundariesFacesIt;
                ++boundariesFacesEntitiesTagsIt;
                ++boundaryConditionIt;
            }
        }


        size_t createVolumeFacesEntities(const int volumeTag,
                                         const int boundariesTags[],
                                         const size_t nBoundaries,
                                         const BoundaryCondition boundariesConditions[],
                                         int& interiorFacesEntity,
                                         size_t& nInteriorFaces,
                                         int boundariesFacesEntytiesTags[],
                                         size_t nBoundariesFaces[],
                                         size_t(*&interiorFacesTetrahedronsIndexes)[2],
                                         uint8_t(*&interiorFacesLocalIndexes)[2],
                                         size_t* boundariesFacesTetrahedronsIndexes[],
                                         uint8_t* boundariesFacesLocalIndexes[])
        {
            size_t* tetrahedronsFacesNodesTagsSeq;
            size_t nTetrahedronsFacesNodesTags;
            GMSHProxy::model::mesh::getFacesByElements(GMSHProxy::model::mesh::TETRAHEDRON_TYPE,
                GMSHProxy::model::mesh::TRIANGLE_FACE_TYPE,
                &tetrahedronsFacesNodesTagsSeq,
                &nTetrahedronsFacesNodesTags,
                volumeTag);

            size_t* tetrahedronsFacesTags;
            GMSHProxy::model::mesh::getFaces(GMSHProxy::model::mesh::TRIANGLE_FACE_TYPE, tetrahedronsFacesNodesTagsSeq, nTetrahedronsFacesNodesTags, &tetrahedronsFacesTags);

            size_t** boundariesFacesTags = new size_t * [nBoundaries];

            getBoundariesFacesTags(boundariesTags, nBoundaries, boundariesFacesTags, nBoundariesFaces);
            allocateMemory(nBoundariesFaces, nBoundaries, boundariesFacesTetrahedronsIndexes, boundariesFacesLocalIndexes);

            size_t nBoundariesTagsSum = std::accumulate(nBoundariesFaces, nBoundariesFaces + nBoundaries, size_t(0));

            size_t nInteriorFaces2x = nTetrahedronsFacesNodesTags / constants::triangle::N_NODES - nBoundariesTagsSum;
            nInteriorFaces = nInteriorFaces2x / 2;
            size_t nFaces = nInteriorFaces + nBoundariesTagsSum;

            interiorFacesTetrahedronsIndexes = new size_t[nInteriorFaces][2];
            interiorFacesLocalIndexes = new uint8_t[nInteriorFaces][2];

            size_t nTetrahedrons = nTetrahedronsFacesNodesTags / (constants::tetrahedron::N_FACES * constants::triangle::N_NODES);
            determineFaceTetrahedronConnections(tetrahedronsFacesTags,
                nTetrahedrons,
                boundariesFacesTags,
                nBoundariesFaces,
                nBoundaries,
                nInteriorFaces,
                nFaces,
                interiorFacesTetrahedronsIndexes,
                interiorFacesLocalIndexes,
                boundariesFacesTetrahedronsIndexes,
                boundariesFacesLocalIndexes);


            GMSHProxy::free(tetrahedronsFacesTags);
            deleteBoundariesFacesTags(nBoundaries, boundariesFacesTags);

            const size_t(*tetrahedronsFacesNodesTags)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES] = (size_t(*)[constants::tetrahedron::N_FACES][constants::triangle::N_NODES])tetrahedronsFacesNodesTagsSeq;
            interiorFacesEntity = createInteriorFacesEntity(tetrahedronsFacesNodesTags, interiorFacesTetrahedronsIndexes, interiorFacesLocalIndexes, nInteriorFaces);

            createBoundariesFacesEntities(tetrahedronsFacesNodesTags,
                boundariesFacesTetrahedronsIndexes,
                boundariesFacesLocalIndexes,
                nBoundariesFaces,
                nBoundaries,
                boundariesConditions,
                boundariesFacesEntytiesTags);

            delete[] boundariesFacesTags;

            GMSHProxy::free(tetrahedronsFacesNodesTagsSeq);

            return nFaces;
        }

        static size_t countFacesEntitiesForRemoving(const int* entitiesTagIt, const size_t nEntities)
        {
            size_t nEntitiesForRemoving = 0;
            for (size_t i = 0; i < nEntities; ++i)
            {
                if (*entitiesTagIt != 0)
                {
                    ++nEntitiesForRemoving;
                }
                ++entitiesTagIt;
            }

            return nEntitiesForRemoving;
        }

        static void fill2DDimTagArray(const int* tagIt,const size_t nTags, int* dimTagIt)
        {
            for (size_t i = 0; i < nTags; ++i)
            {
                if (*tagIt != 0)
                {
                    *dimTagIt = constants::DIMENSION_2D;
                    ++dimTagIt;
                    *dimTagIt = *tagIt;
                    ++dimTagIt;
                }
                ++tagIt;
            }
        }

        void removeCreatedVolumeFacesEntities(const int* entitiesTagIt, const size_t nEntities)
        {
            size_t nEntitiesForRemoving = countFacesEntitiesForRemoving(entitiesTagIt, nEntities);
            size_t nEntitiesForRemovingDimTag = 2 * nEntitiesForRemoving;
            int* removingEntitiesDimTags = new int[nEntitiesForRemovingDimTag];
            fill2DDimTagArray(entitiesTagIt, nEntities, removingEntitiesDimTags);
            
            int ierr;
            gmshModelRemoveEntities(removingEntitiesDimTags, nEntitiesForRemovingDimTag, GMSHProxy::FALSE_C, &ierr);
            if (ierr)
            {
                throwLastError();
            }

            delete[] removingEntitiesDimTags;
        }

        void determineBoundaryConditions(const int* boundaryTagIt, const size_t nBoundaries, BoundaryCondition* conditionTypeIt)
        {
            int ierr;
            for (size_t i = 0; i < nBoundaries; ++i)
            {
                int* boundaryPhysicalTags;
                size_t nBoundaryPhysicalTags;
                gmshModelGetPhysicalGroupsForEntity(constants::DIMENSION_2D, *boundaryTagIt, &boundaryPhysicalTags, &nBoundaryPhysicalTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                if (nBoundaryPhysicalTags == 0)
                {
                    *conditionTypeIt = BoundaryCondition::HOMOGENEOUS_NEWMAN;

                }
                else
                {
                    char* conditionName;
                    gmshModelGetPhysicalName(constants::DIMENSION_2D, *boundaryPhysicalTags, &conditionName, &ierr);
                    if (ierr)
                    {
                        throwLastError();
                    }

                    if (conditionName == nullptr)
                    {
                        *conditionTypeIt = BoundaryCondition::UNDEFINED;
                    }
                    else
                    {
                        switch (*conditionName)
                        {
                        case 'D':
                        {
                            *conditionTypeIt = BoundaryCondition::DIRICHLET;
                            break;
                        }

                        case 'N':
                        {
                            *conditionTypeIt = BoundaryCondition::NEWMAN;
                            break;
                        }

                        case 'C':
                        {
                            *conditionTypeIt = BoundaryCondition::NONCONFORM_INTERAFACE;
                            break;
                        }
                        default:
                        {
                            *conditionTypeIt = BoundaryCondition::UNDEFINED;
                            break;
                        }
                        }
                    }
                    GMSHProxy::free(conditionName);
                }
                GMSHProxy::free(boundaryPhysicalTags);

                ++boundaryTagIt;
                ++conditionTypeIt;
            }
        }

        void determineBoundaryConditions(const int* boundaryTagIt, const size_t nBoundaries, int* conditionTagIt, char** conditionNameIt, BoundaryCondition* conditionTypeIt)
        {
            int ierr;
            for (size_t i = 0; i < nBoundaries; ++i)
            {
                int* boundaryPhysicalTags;
                size_t nBoundaryPhysicalTags;
                gmshModelGetPhysicalGroupsForEntity(constants::DIMENSION_2D, *boundaryTagIt, &boundaryPhysicalTags, &nBoundaryPhysicalTags, &ierr);
                if (ierr)
                {
                    throwLastError();
                }

                if (nBoundaryPhysicalTags == 0)
                {
                    *conditionTagIt = 0;
                    *conditionNameIt = nullptr;
                    *conditionTypeIt = BoundaryCondition::HOMOGENEOUS_NEWMAN;

                }
                else
                {
                    *conditionTagIt = *boundaryPhysicalTags;
                    gmshModelGetPhysicalName(constants::DIMENSION_2D, *conditionTagIt, conditionNameIt, &ierr);
                    if (ierr)
                    {
                        throwLastError();
                    }

                    if (*conditionNameIt == nullptr)
                    {
                        *conditionTypeIt = BoundaryCondition::UNDEFINED;
                    }
                    else
                    {
                        switch (**conditionNameIt)
                        {
                        case 'D':
                        {
                            *conditionTypeIt = BoundaryCondition::DIRICHLET;
                            break;
                        }

                        case 'N':
                        {
                            *conditionTypeIt = BoundaryCondition::NEWMAN;
                            break;
                        }

                        case 'C':
                        {
                            *conditionTypeIt = BoundaryCondition::NONCONFORM_INTERAFACE;
                            break;
                        }
                        default:
                        {
                            *conditionTypeIt = BoundaryCondition::UNDEFINED;
                            break;
                        }
                        }
                    }
                }

                GMSHProxy::free(boundaryPhysicalTags);

                ++boundaryTagIt;
                ++conditionTagIt;
                ++conditionNameIt;
                ++conditionTypeIt;
            }
        }

        void computeInterfacesThermalConductivities(const MaterialPhase* volumeMaterialPhasePtrsIt[],
                                                    const size_t(*intefacesVolumesIndexesIt)[2],
                                                    const size_t nInterfaces,
                                                    double* thermalConductivityIt)
        {
            for (size_t i = 0; i < nInterfaces; ++i)
            {
                *thermalConductivityIt = volumeMaterialPhasePtrsIt[**intefacesVolumesIndexesIt]->thermalConductivity;

                ++intefacesVolumesIndexesIt;
                ++thermalConductivityIt;
            }
        }

        void get3DElementsByPoints(const Coordinates* pointIt, const uint8_t nPoints, size_t* elementTagIt, LocalCoordinates3D* localCoordinatesIt)
        {
            int elementType;
            size_t* nodesTags;
            size_t nNodes;
            int ierr;

            for (uint8_t i = 0; i < nPoints; ++i)
            {
                gmshModelMeshGetElementByCoordinates(pointIt->x, pointIt->y, pointIt->z,
                                                     elementTagIt,
                                                     &elementType,
                                                     &nodesTags,
                                                     &nNodes,
                                                     &(localCoordinatesIt->u), &(localCoordinatesIt->v), &(localCoordinatesIt->w),
                                                     constants::DIMENSION_3D,
                                                     GMSHProxy::FALSE_C,
                                                     &ierr);

                if (ierr)
                {
                    throwLastError();
                }

                ++pointIt;
                ++elementTagIt;
                ++localCoordinatesIt;
            }
        }

        void getModelNodesCoordinates(const int* nodeTagIt, size_t nNodes, Coordinates* coordinatesIt)
        {
            int ierr;
            double xBuf, yBuf, zBuf;
            for (size_t i = 0; i < nNodes; ++i)
            {
                gmshModelGetBoundingBox(constants::DIMENSION_0D, *nodeTagIt, &(coordinatesIt->x), &(coordinatesIt->y), &(coordinatesIt->z), &xBuf, &yBuf, &zBuf, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
                ++nodeTagIt;
                ++coordinatesIt;
            }
        }

        void computeLocalJacobianMatrix(const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT], const double determinant, double localJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT])
        {
            *localJacobianMatrix = (transpJacobianMatrix[4] * transpJacobianMatrix[8] - transpJacobianMatrix[5] * transpJacobianMatrix[7]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[5] * transpJacobianMatrix[6] - transpJacobianMatrix[3] * transpJacobianMatrix[8]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[3] * transpJacobianMatrix[7] - transpJacobianMatrix[4] * transpJacobianMatrix[6]) / determinant;
            ++localJacobianMatrix;

            *localJacobianMatrix = (transpJacobianMatrix[2] * transpJacobianMatrix[7] - transpJacobianMatrix[1] * transpJacobianMatrix[8]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[0] * transpJacobianMatrix[8] - transpJacobianMatrix[2] * transpJacobianMatrix[6]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[1] * transpJacobianMatrix[6] - transpJacobianMatrix[0] * transpJacobianMatrix[7]) / determinant;
            ++localJacobianMatrix;

            *localJacobianMatrix = (transpJacobianMatrix[1] * transpJacobianMatrix[5] - transpJacobianMatrix[2] * transpJacobianMatrix[4]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[2] * transpJacobianMatrix[3] - transpJacobianMatrix[0] * transpJacobianMatrix[5]) / determinant;
            ++localJacobianMatrix;
            *localJacobianMatrix = (transpJacobianMatrix[0] * transpJacobianMatrix[4] - transpJacobianMatrix[1] * transpJacobianMatrix[3]) / determinant;
            ++localJacobianMatrix;
        }

        uint8_t getClosestNodeIndex(const Coordinates& node, const size_t* nodeTagIt, const uint8_t nNodes)
        {
            double* nodeCoordinates;
            size_t nCoordinates;
            double* parametricCoordinates;
            size_t nParametricCoordinates;

            int dim, tag, ierr;

            uint8_t index = 0;
            double errMax;

            if (nNodes != 0)
            {
                gmshModelMeshGetNode(*nodeTagIt, &nodeCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, &dim, &tag, &ierr);
                if (ierr)
                {
                    throwLastError();
                }
                errMax = abs(node.x - nodeCoordinates[0]) + abs(node.y - nodeCoordinates[1]) + abs(node.z - nodeCoordinates[2]);

                gmshFree(nodeCoordinates);
                gmshFree(parametricCoordinates);

                ++nodeTagIt;
                index = 0;
            }
            else
            {
                return nNodes;
            }

            for (uint8_t i = 1; i < nNodes; ++i)
            {
                gmshModelMeshGetNode(*nodeTagIt, &nodeCoordinates, &nCoordinates, &parametricCoordinates, &nParametricCoordinates, &dim, &tag, &ierr);

                if (ierr)
                {
                    throwLastError();
                }

                double err = abs(node.x - nodeCoordinates[0]) + abs(node.y - nodeCoordinates[1]) + abs(node.z - nodeCoordinates[2]);

                if (err < errMax)
                {
                    index = i;
                    errMax = err;
                }
                gmshFree(nodeCoordinates);
                gmshFree(parametricCoordinates);

                ++nodeTagIt;
            }
        }

        void translateTo0FaceLocalPoints(const double* integrationLocalCoordinateIt, const size_t nSteps, LocalCoordinates3D* localPointIt)
        {
            for (size_t i = 0; i < nSteps; ++i)
            {
                localPointIt->v = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->u = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->w = 0.0;
                ++integrationLocalCoordinateIt;

                ++localPointIt;
            }
        }

        void translateTo1FaceLocalPoints(const double* integrationLocalCoordinateIt, const size_t nSteps, LocalCoordinates3D* localPointIt)
        {
            for (size_t i = 0; i < nSteps; ++i)
            {
                localPointIt->u = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->w = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->v = 0.0;
                ++integrationLocalCoordinateIt;

                ++localPointIt;
            }
        }

        void translateTo2FaceLocalPoints(const double* integrationLocalCoordinateIt, const size_t nSteps, LocalCoordinates3D* localPointIt)
        {
            for (size_t i = 0; i < nSteps; ++i)
            {
                localPointIt->w = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->v = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->u = 0.0;
                ++integrationLocalCoordinateIt;

                ++localPointIt;
            }
        }

        void translateTo3FaceLocalPoints(const double* integrationLocalCoordinateIt, const size_t nSteps, LocalCoordinates3D* localPointIt)
        {
            for (size_t i = 0; i < nSteps; ++i)
            {
                localPointIt->u = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->v = *integrationLocalCoordinateIt;
                ++integrationLocalCoordinateIt;
                localPointIt->w = 1.0 - localPointIt->u - localPointIt->v;
                ++integrationLocalCoordinateIt;

                ++localPointIt;
            }
        }
        


        void getIntegrationByTetrahedronFace(const char methodName[], LocalCoordinates3D* facesPoints[constants::tetrahedron::N_FACES], double* &weights, size_t &nSteps)
        {
            int ierr;
            double* localCoordinates;
            size_t nLocalCoords;
            size_t nIntegrationStepsProxy;
            gmshModelMeshGetIntegrationPoints(GMSHProxy::model::mesh::TRIANGLE_TYPE, methodName, &localCoordinates, &nLocalCoords, &weights, &nSteps, &ierr);

            if (ierr)
            {
                throwLastError();
            }

            facesPoints[0] = new LocalCoordinates3D[constants::tetrahedron::N_FACES * nSteps];
            facesPoints[1] = facesPoints[0] + nSteps;
            facesPoints[2] = facesPoints[1] + nSteps;
            facesPoints[3] = facesPoints[2] + nSteps;

            translateTo0FaceLocalPoints(localCoordinates, nSteps, facesPoints[0]);
            translateTo1FaceLocalPoints(localCoordinates, nSteps, facesPoints[1]);
            translateTo2FaceLocalPoints(localCoordinates, nSteps, facesPoints[2]);
            translateTo3FaceLocalPoints(localCoordinates, nSteps, facesPoints[3]);

            gmshFree(localCoordinates);

        }

        void freeIntegration(LocalCoordinates3D* facesPoints[constants::tetrahedron::N_FACES], double* weights)
        {
            delete[] *facesPoints;
            gmshFree(weights);
        }
};
