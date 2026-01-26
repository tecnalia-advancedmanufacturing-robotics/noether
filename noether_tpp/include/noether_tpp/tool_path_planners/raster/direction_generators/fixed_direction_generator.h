#pragma once

#include <noether_tpp/tool_path_planners/raster/raster_planner.h>
#include <noether_tpp/tool_path_planners/cross_raster/cross_raster_planner.h>
#include <noether_tpp/tool_path_planners/grid_raster/grid_raster_planner.h>

namespace noether
{
/**
 * @ingroup direction_generators
 * @brief Generates the raster direction based on a fixed input
 */
class FixedDirectionGenerator : public DirectionGenerator, public CrossDirectionGenerator, public GridDirectionGenerator
{
public:
  FixedDirectionGenerator(const Eigen::Vector3d& direction);
  Eigen::Vector3d generate(const pcl::PolygonMesh& mesh) const override final;

private:
  const Eigen::Vector3d direction_;
};

}  // namespace noether
