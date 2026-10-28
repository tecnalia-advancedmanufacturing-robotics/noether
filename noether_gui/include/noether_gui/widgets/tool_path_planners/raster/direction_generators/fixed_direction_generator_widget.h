#pragma once

#include <noether_gui/widgets.h>

#include <noether_tpp/tool_path_planners/raster/raster_planner.h>

namespace Ui
{
class Vector3dEditor;
}

namespace noether
{
class FixedDirectionGeneratorWidget : public DirectionGeneratorWidget
{
public:
  FixedDirectionGeneratorWidget(QWidget* parent = nullptr);

  DirectionGenerator::ConstPtr create() const override;

  void configure(const YAML::Node&) override;
  void save(YAML::Node&) const override;

private:
  Ui::Vector3dEditor* cut_direction_ui_;
  Ui::Vector3dEditor* cut_normal_ui_;
};

}  // namespace noether
