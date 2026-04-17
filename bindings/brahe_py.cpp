#include "brahe/body_system.h"
#include "brahe/conics.h"
#include "brahe/event_detector.h"
#include "brahe/patcher.h"
#include "brahe/trajectory.h"
#include "brahe/trajectory_builder.h"
#include "brahe/two_body.h"
#include "brahe/types.h"

#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace brahe {
namespace {

class PyEventDetector {
public:
    PyEventDetector(std::shared_ptr<BodySystem> bodies, Tolerances tolerances)
        : bodies_(std::move(bodies)), tolerances_(tolerances) {}

    py::object find_next_event(const EventSearchRequest& req, bool diagnostics) const {
        PredictedEvent event{};
        if (diagnostics) {
            EventDetectorDiagnostics diag{};
            EventDetector detector(*bodies_, tolerances_, &diag);
            SolveStatus status = detector.find_next_event(req, event);
            return py::make_tuple(status, event, diag);
        }

        EventDetector detector(*bodies_, tolerances_);
        SolveStatus status = detector.find_next_event(req, event);
        return py::make_tuple(status, event);
    }

private:
    std::shared_ptr<BodySystem> bodies_;
    Tolerances tolerances_;
};

class PyTrajectoryBuilder {
public:
    PyTrajectoryBuilder(std::shared_ptr<BodySystem> bodies, Tolerances tolerances)
        : bodies_(std::move(bodies)), tolerances_(tolerances) {}

    py::tuple build_preview(const PreviewRequest& req) const {
        Trajectory trajectory{};
        TrajectoryBuilder builder(*bodies_, tolerances_);
        SolveStatus status = builder.build_preview(req, trajectory);
        return py::make_tuple(status, trajectory);
    }

private:
    std::shared_ptr<BodySystem> bodies_;
    Tolerances tolerances_;
};

class PyPatcher {
public:
    PyPatcher(std::shared_ptr<BodySystem> bodies, Tolerances tolerances)
        : bodies_(std::move(bodies)), tolerances_(tolerances) {}

    py::tuple patch_parent_to_child(double event_time, BodyId parent_body,
                                    BodyId child_body,
                                    const State2& spacecraft_in_parent_frame) const {
        State2 out{};
        Patcher patcher(*bodies_, tolerances_);
        SolveStatus status = patcher.patch_parent_to_child(
            event_time, parent_body, child_body, spacecraft_in_parent_frame, out);
        return py::make_tuple(status, out);
    }

    py::tuple patch_child_to_parent(double event_time, BodyId child_body,
                                    const State2& spacecraft_in_child_frame) const {
        State2 out{};
        Patcher patcher(*bodies_, tolerances_);
        SolveStatus status = patcher.patch_child_to_parent(
            event_time, child_body, spacecraft_in_child_frame, out);
        return py::make_tuple(status, out);
    }

private:
    std::shared_ptr<BodySystem> bodies_;
    Tolerances tolerances_;
};

void append_body_ids(const BodySystem& bodies, BodyId body, std::vector<BodyId>& ids) {
    ids.push_back(body);
    const BodyId* begin = bodies.children_begin(body);
    const BodyId* end = bodies.children_end(body);
    for (const BodyId* p = begin; p != end; ++p) {
        append_body_ids(bodies, *p, ids);
    }
}

}  // namespace
}  // namespace brahe

PYBIND11_MODULE(_brahe, m) {
    using namespace brahe;

    m.doc() = "Python bindings for Brahe patched-conic trajectory tools";
    m.attr("InvalidBody") = py::int_(InvalidBody);

    py::enum_<ConicType>(m, "ConicType")
        .value("Ellipse", ConicType::Ellipse)
        .value("ParabolaLike", ConicType::ParabolaLike)
        .value("Hyperbola", ConicType::Hyperbola);

    py::enum_<EventType>(m, "EventType")
        .value("None_", EventType::None)
        .value("Impact", EventType::Impact)
        .value("SoiEntry", EventType::SoiEntry)
        .value("SoiExit", EventType::SoiExit)
        .value("Burn", EventType::Burn)
        .value("TimeLimit", EventType::TimeLimit);

    py::enum_<SolveStatus>(m, "SolveStatus")
        .value("Ok", SolveStatus::Ok)
        .value("InvalidInput", SolveStatus::InvalidInput)
        .value("NoConvergence", SolveStatus::NoConvergence)
        .value("NoSolution", SolveStatus::NoSolution)
        .value("NumericalFailure", SolveStatus::NumericalFailure)
        .value("CapacityExceeded", SolveStatus::CapacityExceeded);

    py::enum_<LambertDirection>(m, "LambertDirection")
        .value("CCW", LambertDirection::CCW)
        .value("CW", LambertDirection::CW);

    py::enum_<BurnFrame>(m, "BurnFrame")
        .value("Inertial", BurnFrame::Inertial)
        .value("ProgradeRadial", BurnFrame::ProgradeRadial)
        .value("ProgradeNormal", BurnFrame::ProgradeNormal);

    py::class_<Vec2>(m, "Vec2")
        .def(py::init<double, double>(), py::arg("x") = 0.0, py::arg("y") = 0.0)
        .def_readwrite("x", &Vec2::x)
        .def_readwrite("y", &Vec2::y)
        .def("__repr__", [](const Vec2& v) {
            return "Vec2(x=" + std::to_string(v.x) + ", y=" + std::to_string(v.y) + ")";
        });

    py::class_<State2>(m, "State2")
        .def(py::init<Vec2, Vec2>(), py::arg("r") = Vec2{}, py::arg("v") = Vec2{})
        .def_readwrite("r", &State2::r)
        .def_readwrite("v", &State2::v);

    py::class_<Tolerances>(m, "Tolerances")
        .def(py::init<>())
        .def_readwrite("position_epsilon", &Tolerances::position_epsilon)
        .def_readwrite("velocity_epsilon", &Tolerances::velocity_epsilon)
        .def_readwrite("angle_epsilon", &Tolerances::angle_epsilon)
        .def_readwrite("time_epsilon", &Tolerances::time_epsilon)
        .def_readwrite("root_epsilon", &Tolerances::root_epsilon)
        .def_readwrite("lambert_residual_epsilon", &Tolerances::lambert_residual_epsilon)
        .def_readwrite("parabolic_eccentricity_band",
                       &Tolerances::parabolic_eccentricity_band)
        .def_readwrite("max_event_refine_iterations",
                       &Tolerances::max_event_refine_iterations);

    py::class_<BodyDef>(m, "BodyDef")
        .def(py::init<>())
        .def_readwrite("id", &BodyDef::id)
        .def_readwrite("parent_id", &BodyDef::parent_id)
        .def_readwrite("mu", &BodyDef::mu)
        .def_readwrite("radius", &BodyDef::radius)
        .def_readwrite("soi_radius", &BodyDef::soi_radius)
        .def_readwrite("orbit_radius", &BodyDef::orbit_radius)
        .def_readwrite("angular_rate", &BodyDef::angular_rate)
        .def_readwrite("phase_at_epoch", &BodyDef::phase_at_epoch);

    py::class_<ConicElements2D>(m, "ConicElements2D")
        .def(py::init<>())
        .def_readwrite("type", &ConicElements2D::type)
        .def_readwrite("mu", &ConicElements2D::mu)
        .def_readwrite("semi_major_axis", &ConicElements2D::semi_major_axis)
        .def_readwrite("eccentricity", &ConicElements2D::eccentricity)
        .def_readwrite("semi_latus_rectum", &ConicElements2D::semi_latus_rectum)
        .def_readwrite("angular_momentum_z", &ConicElements2D::angular_momentum_z)
        .def_readwrite("argument_of_periapsis", &ConicElements2D::argument_of_periapsis)
        .def_readwrite("true_anomaly", &ConicElements2D::true_anomaly)
        .def_readwrite("mean_anomaly", &ConicElements2D::mean_anomaly)
        .def_readwrite("eccentric_anomaly", &ConicElements2D::eccentric_anomaly)
        .def_readwrite("hyperbolic_anomaly", &ConicElements2D::hyperbolic_anomaly);

    py::class_<Segment>(m, "Segment")
        .def(py::init<>())
        .def_readwrite("central_body", &Segment::central_body)
        .def_readwrite("start_time", &Segment::start_time)
        .def_readwrite("initial_state", &Segment::initial_state)
        .def_readwrite("end_time", &Segment::end_time)
        .def_readwrite("end_reason", &Segment::end_reason);

    py::class_<Trajectory>(m, "Trajectory")
        .def(py::init<>())
        .def_readwrite("segments", &Trajectory::segments);

    py::class_<PreviewRequest>(m, "PreviewRequest")
        .def(py::init<>())
        .def_readwrite("central_body", &PreviewRequest::central_body)
        .def_readwrite("start_time", &PreviewRequest::start_time)
        .def_readwrite("initial_state", &PreviewRequest::initial_state)
        .def_readwrite("end_time", &PreviewRequest::end_time)
        .def_readwrite("max_segments", &PreviewRequest::max_segments);

    py::class_<PredictedEvent>(m, "PredictedEvent")
        .def(py::init<>())
        .def_readwrite("type", &PredictedEvent::type)
        .def_readwrite("time", &PredictedEvent::time)
        .def_readwrite("from_body", &PredictedEvent::from_body)
        .def_readwrite("to_body", &PredictedEvent::to_body)
        .def_readwrite("state", &PredictedEvent::state);

    py::class_<EventSearchRequest>(m, "EventSearchRequest")
        .def(py::init<>())
        .def_readwrite("central_body", &EventSearchRequest::central_body)
        .def_readwrite("start_time", &EventSearchRequest::start_time)
        .def_readwrite("initial_state", &EventSearchRequest::initial_state)
        .def_readwrite("time_limit", &EventSearchRequest::time_limit);

    py::class_<EventDetectorDiagnostics>(m, "EventDetectorDiagnostics")
        .def(py::init<>())
        .def_readwrite("scan_steps", &EventDetectorDiagnostics::scan_steps)
        .def_readwrite("root_refinements", &EventDetectorDiagnostics::root_refinements)
        .def_readwrite("impact_minimum_refinements",
                       &EventDetectorDiagnostics::impact_minimum_refinements)
        .def_readwrite("soi_exit_minimum_refinements",
                       &EventDetectorDiagnostics::soi_exit_minimum_refinements)
        .def_readwrite("child_entry_minimum_refinements",
                       &EventDetectorDiagnostics::child_entry_minimum_refinements);

    py::class_<BodySystem, std::shared_ptr<BodySystem>>(m, "BodySystem")
        .def("get_body", [](const BodySystem& self, BodyId id) -> py::object {
            const BodyDef* body = self.get_body(id);
            if (body == nullptr) return py::none();
            return py::cast(*body);
        })
        .def("position_in_parent", &BodySystem::position_in_parent)
        .def("velocity_in_parent", &BodySystem::velocity_in_parent)
        .def("state_in_ancestor_frame",
             [](const BodySystem& self, BodyId body, BodyId ancestor, double t) {
                 State2 out{};
                 SolveStatus status = self.state_in_ancestor_frame(body, ancestor, t, out);
                 return py::make_tuple(status, out);
             })
        .def("state_in_root_frame", &BodySystem::state_in_root_frame)
        .def("root_id", &BodySystem::root_id)
        .def("body_count", &BodySystem::body_count)
        .def("depth", &BodySystem::depth)
        .def("children", [](const BodySystem& self, BodyId id) {
            std::vector<BodyId> children;
            const BodyId* begin = self.children_begin(id);
            const BodyId* end = self.children_end(id);
            children.assign(begin, end);
            return children;
        })
        .def("body_ids", [](const BodySystem& self) {
            std::vector<BodyId> ids;
            BodyId root = self.root_id();
            if (root != InvalidBody) append_body_ids(self, root, ids);
            return ids;
        });

    py::class_<BodySystemBuilder>(m, "BodySystemBuilder")
        .def(py::init<>())
        .def("add_body", &BodySystemBuilder::add_body)
        .def("build", [](const BodySystemBuilder& self) -> py::tuple {
            BodySystem out{};
            SolveStatus status = self.build(out);
            if (status != SolveStatus::Ok) {
                return py::make_tuple(status, py::none());
            }
            return py::make_tuple(status, py::cast(std::make_shared<BodySystem>(out)));
        });

    py::class_<TwoBody>(m, "TwoBody")
        .def_static("classify", [](double mu, const State2& state) {
            ConicType type{};
            SolveStatus status = TwoBody::classify(mu, state, type);
            return py::make_tuple(status, type);
        })
        .def_static("to_elements", [](double mu, const State2& state) {
            ConicElements2D elements{};
            SolveStatus status = TwoBody::to_elements(mu, state, elements);
            return py::make_tuple(status, elements);
        })
        .def_static("from_elements", [](const ConicElements2D& elements) {
            State2 state{};
            SolveStatus status = TwoBody::from_elements(elements, state);
            return py::make_tuple(status, state);
        })
        .def_static("propagate", [](double mu, const State2& initial, double dt) {
            State2 state{};
            SolveStatus status = TwoBody::propagate(mu, initial, dt, state);
            return py::make_tuple(status, state);
        })
        .def_static("specific_energy", &TwoBody::specific_energy)
        .def_static("specific_angular_momentum_z",
                    &TwoBody::specific_angular_momentum_z)
        .def_static("eccentricity", &TwoBody::eccentricity)
        .def_static("semi_major_axis", &TwoBody::semi_major_axis)
        .def_static("semi_latus_rectum", &TwoBody::semi_latus_rectum)
        .def_static("periapsis_radius", &TwoBody::periapsis_radius)
        .def_static("apoapsis_radius", &TwoBody::apoapsis_radius);

    py::class_<PyEventDetector>(m, "EventDetector")
        .def(py::init<std::shared_ptr<BodySystem>, Tolerances>(),
             py::arg("bodies"), py::arg("tolerances") = Tolerances{})
        .def("find_next_event", &PyEventDetector::find_next_event,
             py::arg("request"), py::arg("diagnostics") = false);

    py::class_<PyTrajectoryBuilder>(m, "TrajectoryBuilder")
        .def(py::init<std::shared_ptr<BodySystem>, Tolerances>(),
             py::arg("bodies"), py::arg("tolerances") = Tolerances{})
        .def("build_preview", &PyTrajectoryBuilder::build_preview);

    py::class_<PyPatcher>(m, "Patcher")
        .def(py::init<std::shared_ptr<BodySystem>, Tolerances>(),
             py::arg("bodies"), py::arg("tolerances") = Tolerances{})
        .def("patch_parent_to_child", &PyPatcher::patch_parent_to_child)
        .def("patch_child_to_parent", &PyPatcher::patch_child_to_parent);
}
