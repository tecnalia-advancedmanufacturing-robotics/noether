#pragma once

#include <noether_tpp/tool_path_planners/raster/raster_planner.h>
#include <noether_tpp/tool_path_planners/cross_raster/cross_raster_planner.h>

namespace noether
{
/**
 * @ingroup direction_generators
 * @brief Generates the raster direction based on a fixed input
 */
class FixedDirectionGenerator : public DirectionGenerator, public CrossDirectionGenerator
{
public:
  FixedDirectionGenerator(const Eigen::Vector3d& cut_direction, const Eigen::Vector3d& cut_normal);
  Eigen::Vector3d generate_cut_direction(const pcl::PolygonMesh& mesh) const override final;
  Eigen::Vector3d generate_mesh_normal(const pcl::PolygonMesh& mesh) const override final;

private:
  const Eigen::Vector3d cut_direction_;
  const Eigen::Vector3d cut_normal_;
};

}  // namespace noether
