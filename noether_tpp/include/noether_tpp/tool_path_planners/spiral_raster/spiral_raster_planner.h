/**
 * @file spiral_raster_planner.h
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

#include <functional>
#include <memory>

#include <noether_tpp/core/tool_path_planner.h>

namespace noether
{
/**
 * @brief Supplies the in-plane x axis and the surface normal used to lay out a spiral.
 */
struct SpiralDirectionGenerator
{
  using Ptr = std::unique_ptr<SpiralDirectionGenerator>;
  using ConstPtr = std::unique_ptr<const SpiralDirectionGenerator>;

  virtual ~SpiralDirectionGenerator() = default;
  virtual Eigen::Vector3d generate_cut_direction(const pcl::PolygonMesh& mesh) const = 0;
  virtual Eigen::Vector3d generate_mesh_normal(const pcl::PolygonMesh& mesh) const = 0;
};

/**
 * @brief Plans one spiral on a mesh. Subclasses generate the waypoints; this class aligns their x-axis with travel.
 */
class SpiralRasterPlanner : public ToolPathPlanner
{
public:
  explicit SpiralRasterPlanner(SpiralDirectionGenerator::ConstPtr dir_gen);

  ToolPaths plan(const pcl::PolygonMesh& mesh) const override final;

  void setPointSpacing(double point_spacing);
  void setLineSpacing(double line_spacing);
  void setMinHoleSize(double min_hole_size);

protected:
  virtual ToolPaths planImpl(const pcl::PolygonMesh& mesh) const = 0;

  SpiralDirectionGenerator::ConstPtr dir_gen_;

  /** @brief Distance between waypoints along the spiral (m) */
  double point_spacing_;
  /** @brief Radial distance between successive turns (m) */
  double line_spacing_;
  /** @brief Gap above which a missed projection splits the spiral (m) */
  double min_hole_size_;
};

struct SpiralRasterPlannerFactory : public ToolPathPlannerFactory
{
  double point_spacing;
  double line_spacing;
  double min_hole_size;
  std::function<std::unique_ptr<const SpiralDirectionGenerator>()> direction_gen;
};

}  // namespace noether
