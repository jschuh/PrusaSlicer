#ifndef slic3r_ZDither_hpp
#define slic3r_Zdither_hpp

#include "Point.hpp"
#include "TriangleMeshSlicer.hpp"

namespace Slic3r {

    struct SubLayers
{
    ExPolygons bottom_; // polygons filling the botthom 0.25% of layer thickness
    ExPolygons halfUp_; // polygons filling only half a layer leaving 0.25% at top and botton not filled, located above bottom_
    ExPolygons halfDn_; // similar to halfUp_ but located under top_
    ExPolygons top_; // polygons filling the top 0.25% of layer thickness
    SubLayers()                       = default;
    SubLayers(const SubLayers &other) = default;
    SubLayers(SubLayers &&other)      = default;
    bool empty() const { return bottom_.empty() && halfUp_.empty() && halfDn_.empty() && top_.empty(); };
};

void midslice_zs(const indexed_triangle_set &mesh,
                 const std::vector<float> &  zs,
                 const Transform3d &         trafo,
                 float                       nozzle_diameter,
                 std::vector<float> *        mid_zs,
                 std::vector<int> *          upwrd_mididx,
                 std::vector<int> *          dnwrd_mididx);

std::vector<ExPolygons> apply_z_dither(std::vector<ExPolygons> &layers,
                                       std::vector<ExPolygons> &mid_layers,
                                       const std::vector<bool> &do_low,
                                       const std::vector<bool> &do_high,
                                       std::vector<SubLayers> * sublayers);

std::vector<ExPolygons> z_dither(const indexed_triangle_set &mesh,
                                 const std::vector<float> &  zs,
                                 const MeshSlicingParamsEx &params,
                                 std::vector<ExPolygons> &    layers,
                                 std::vector<SubLayers> *     sublayers,
                                 const std::function<void()> &throw_on_cancel_callback);

} // namespace Slic3r
#endif
