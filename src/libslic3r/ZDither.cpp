#include "ZDither.hpp"
#include "ExPolygon.hpp"
#include "TriangleMeshSlicer.hpp"


namespace Slic3r {

// Find z where additional midlevel slices need to be made for z-dithering
// It is assumed that original layers corresponding to zs[i] will be filled
// from zs[i] - layerHight/2 till zs[i] + layerHeight/2
// Required midlevel slices are to be computed halfway between layers,
// starting with top of layer 0 and ending with top of last layer ;

void midslice_zs(const indexed_triangle_set &mesh,
                 const std::vector<float> &  zs,
                 const Transform3d &         trafo,
                 float                       nozzle_diameter,
                 std::vector<float> *        mid_zs,
                 std::vector<int> *         upwrd_mididx,
                 std::vector<int> *         dnwrd_mididx)
{
    mid_zs->clear();
    upwrd_mididx->clear();
    dnwrd_mididx->clear();

    auto n_zs = zs.size();
    if (n_zs < 2) return;

    std::vector<stl_vertex> transformed_vertices(mesh.vertices);
    auto                    trafo_float = trafo.cast<float>();
    for (stl_vertex &v : transformed_vertices) { v = trafo_float * v; }

    std::vector<float> candidate_zs(n_zs);
    std::vector<float> layerHeight(n_zs);
    float              layerHeightMax = 0;

    for (auto i = 0; i < n_zs - 1; i++) {
        // This may be a bit inaccurate if layer heights vary, but will still be an improvement.
        candidate_zs[i] = (zs[i + 1] + zs[i]) / 2;
        layerHeight[i]  = zs[i + 1] - zs[i];
        layerHeightMax  = std::max(layerHeight[i], layerHeightMax);
    }
    layerHeight[n_zs - 1]  = layerHeight[n_zs - 2];
    candidate_zs[n_zs - 1] = zs[n_zs - 1] + layerHeight[n_zs - 1] / 2;

    upwrd_mididx->assign(n_zs, -1);
    dnwrd_mididx->assign(n_zs, -1);

    float nDia = 1.0; // Cooficient to avoid cutting expolys with too vertical triangles (see below)
    // With nozzle_diameter = 0.4 and layer_height = 0.25
    // nDia = 1.5 results in slopes below 22.5 degrees; nDia = 1 results in slopes under 32 degrees
    // nDia = 0.7 - 42 degrees; nDia = 0.5 51.4 degrees

    float eps  = layerHeightMax / 100;

    for (auto i = 0; i < mesh.indices.size(); ++i) {
        const stl_triangle_vertex_indices &facetVertices = mesh.indices[i];
        stl_vertex vertex[3] = {transformed_vertices[facetVertices[0]], transformed_vertices[facetVertices[1]],
                                transformed_vertices[facetVertices[2]]};

        Vec3f norm((vertex[1] - vertex[0]).cross(vertex[2] - vertex[0]).normalized());
        // Skip triangle if over a thickest layer it would not fit nDia * nozzle_diameter
        if (nozzle_diameter * sqrt(1 - norm[2] * norm[2]) * nDia  > fabs(norm[2]) * layerHeightMax ) continue;

         float z_min = std::min(std::min(vertex[0][2], vertex[1][2]), vertex[2][2]);
         float z_max = std::max(std::max(vertex[0][2], vertex[1][2]), vertex[2][2]);

         int start = std::lower_bound(candidate_zs.begin(), candidate_zs.end(), z_min - eps) - candidate_zs.begin();

         for (auto cc = start; cc < candidate_zs.size() && candidate_zs[cc] < z_max + eps; cc++) {
            // Ignore facets that overhang or are too vertical to fit nDia nozzles at layer height
            if (nozzle_diameter * sqrt(1 - norm[2] * norm[2]) * nDia < norm[2] * layerHeight[cc]) {
                if (norm[2] > 0)
                    upwrd_mididx->at(cc) = cc;
                else
                    dnwrd_mididx->at(cc) = cc;
            }
        }
    }

    if (std::all_of(upwrd_mididx->begin(), upwrd_mididx->end(), [](int idx) { return idx == -1; }) &&
        std::all_of(dnwrd_mididx->begin(), dnwrd_mididx->end(), [](int idx) { return idx == -1; }))
        return;

    // Retrn mid_zs which will contribute to z-dithering
    for (int i = 0; i < n_zs; i++) {
        if (upwrd_mididx->at(i) + dnwrd_mididx->at(i) != -2) {
            int count = mid_zs->size();
            mid_zs->push_back(candidate_zs[i]);
            if (upwrd_mididx->at(i) != -1) upwrd_mididx->at(i) = count;
            if (dnwrd_mididx->at(i) != -1) dnwrd_mididx->at(i) = count;
        }
    }
    return;
}

std::vector<ExPolygons> apply_z_dither(std::vector<ExPolygons> &expolys,
                                       std::vector<ExPolygons> &expolys_mid,
                                       const std::vector<int> & upwrd_mididx,
                                       const std::vector<int> & dnwrd_mididx,
                                       std::vector<SubLayers> * sublayers)
{
    sublayers->clear();
    sublayers->resize(expolys.size());
    std::vector<ExPolygons> out(expolys.size());

    out[0] = std::move(expolys[0]);     // Do not make sublayers of first layer
    for (auto ll = 1; ll < expolys.size(); ll++) {
        // idx0 - bottom of layer, idx1 - top of layer
        int  upwrd_idx0  = upwrd_mididx[ll - 1];
        int  dnwrd_idx0  = dnwrd_mididx[ll - 1];
        int  upwrd_idx1  = upwrd_mididx[ll];
        int  dnwrd_idx1  = dnwrd_mididx[ll];

        auto useMidCut = [](int idx) { return idx != -1; };

        if (!useMidCut(upwrd_idx0) && !useMidCut(dnwrd_idx0) && !useMidCut(upwrd_idx1) && !useMidCut(dnwrd_idx1)) {
            out[ll] = std::move(expolys[ll]);
            continue;
        }

        ExPolygons bottom, middleUp, middleDn, top;

        if (useMidCut(upwrd_idx0) || useMidCut(upwrd_idx1)) {
            bottom = std::move(
                diff_ex(useMidCut(upwrd_idx0) ? expolys_mid[upwrd_idx0] : expolys[ll],
                        useMidCut(upwrd_idx1) ? expolys_mid[upwrd_idx1] : expolys[ll]));
        }
        if (useMidCut(upwrd_idx1)) {
            middleUp = std::move(diff_ex(expolys[ll], expolys_mid[upwrd_idx1]));
        }

        if (useMidCut(dnwrd_idx0) || useMidCut(dnwrd_idx1)) {
            top = std::move(
                diff_ex(useMidCut(dnwrd_idx1) ? expolys_mid[dnwrd_idx1] : expolys[ll],
                        useMidCut(dnwrd_idx0) ? expolys_mid[dnwrd_idx0] : expolys[ll]));
        }
        if (useMidCut(dnwrd_idx0)) {
            middleDn = std::move(diff_ex(expolys[ll], expolys_mid[dnwrd_idx0]));
        }

        ExPolygons whole;
        if (useMidCut(upwrd_idx1) && useMidCut(dnwrd_idx0))
            whole = std::move(intersection_ex(expolys_mid[dnwrd_idx0], expolys_mid[upwrd_idx1]));
        else if (useMidCut(upwrd_idx1))
            whole = std::move(intersection_ex(expolys[ll], expolys_mid[upwrd_idx1]));
        else if (useMidCut(dnwrd_idx0))
            whole = std::move(intersection_ex(expolys[ll], expolys_mid[dnwrd_idx0]));
        else {
            out[ll] = std::move(expolys[ll]);
            continue;
        }
        out[ll] = std::move(whole);
        if (bottom.empty() != middleUp.empty() || middleDn.empty() != top.empty()) {
            BOOST_LOG_TRIVIAL(error)
                << "z-dithering: internal error";
        } else {
            sublayers->at(ll).bottom_ = std::move(bottom);
            sublayers->at(ll).halfUp_ = std::move(middleUp);
            sublayers->at(ll).halfDn_ = std::move(middleDn);
            sublayers->at(ll).top_    = std::move(top);
        }
    }

    return out;
}


std::vector<ExPolygons> z_dither(const indexed_triangle_set &mesh,
                                 const std::vector<float> &zs,
                                 const MeshSlicingParamsEx &params,
                                 std::vector<ExPolygons> &    expolys,
                                 std::vector<SubLayers> *     sublayers,
                                 const std::function<void()> &throw_on_cancel_callback)
{
    std::vector<float> mid_zs;
    std::vector<int>   upwrd_mididx;
    std::vector<int>   dnwrd_mididx;
    midslice_zs(mesh, zs, params.trafo, params.nozzle_diameter, &mid_zs, &upwrd_mididx, &dnwrd_mididx);
    if (!mid_zs.empty()) {
        std::vector<ExPolygons> expolys_mid = slice_mesh_ex(mesh, mid_zs, params, throw_on_cancel_callback);
        return apply_z_dither(expolys, expolys_mid, upwrd_mididx, dnwrd_mididx, sublayers);
    }

    *sublayers = std::vector<SubLayers>(expolys.size());
    return expolys;
}

} // namespace Slic3r