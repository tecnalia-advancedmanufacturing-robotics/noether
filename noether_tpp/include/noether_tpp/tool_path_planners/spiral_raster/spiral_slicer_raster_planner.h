/**
 * @file spiral_slicer_raster_planner.h
 * @copyright Copyright (c) 2021, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (Apache License)
 * @par
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * @par
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

#include <noether_tpp/tool_path_planners/spiral_raster/spiral_raster_planner.h>

namespace noether
{
/**
 * @brief Generates one Archimedean spiral on a mesh.
 * @details The spiral is centered on the mesh centroid. delta_x and delta_y shift that center in the plane spanned by
 * the cut direction and the in-plane axis perpendicular to the mesh normal. Each sample is projected onto the mesh
 * along the mesh normal, and the waypoint z-axis is that surface normal. line_spacing is the radial pitch.
 * point_spacing is the distance along the spiral.
 */
class SpiralSlicerRasterPlanner : public SpiralRasterPlanner
{
public:
  explicit SpiralSlicerRasterPlanner(SpiralDirectionGenerator::ConstPtr dir_gen);

  void setSearchRadius(double search_radius);
  void setMinSegmentSize(double min_segment_size);
  void setCenterOffset(double delta_x, double delta_y);

protected:
  ToolPaths planImpl(const pcl::PolygonMesh& mesh) const override;

  /** @brief Search radius for averaging vertex normals (m). Zero uses the intersected face normal. */
  double search_radius_;
  /** @brief Segments shorter than this are dropped (m) */
  double min_segment_size_;
  /** @brief Shift of the spiral center from the centroid along the cut direction (m) */
  double delta_x_;
  /** @brief Shift of the spiral center from the centroid along the in-plane y axis (m) */
  double delta_y_;
};

struct SpiralSlicerRasterPlannerFactory : public SpiralRasterPlannerFactory
{
  double min_segment_size;
  double search_radius;
  double delta_x;
  double delta_y;

  ToolPathPlanner::ConstPtr create() const override;
};

}  // namespace noether
