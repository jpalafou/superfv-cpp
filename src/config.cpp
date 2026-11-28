#include "config.hpp"

#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "io.hpp"

namespace {
    struct Config {
        // required config members start with invalid values.
        std::string         integrator;
        int                 nghost = -1;
        int                 nx = -1;
        int                 p = -1;
        std::string         setup;
        Real                tstop = -1.0;

        // optional config members and their defaults.
        Real                cfl = 0.8;
        bool                clip_last_dt = true;
        std::string         fallback_cascade = "MUSCL0";
        std::size_t         fprint = 0;
        Real                gamma = 1.4;
        std::size_t         max_steps =  std::numeric_limits<std::size_t>::max();
        std::string         MUSCL_limiter = "moncen";
        bool                NAD = true;
        std::size_t         nout = 0;
        std::size_t         nrev = 1;
        std::size_t         ny = 1;
        std::size_t         nz = 1;
        std::string         output_dir = "out";
        bool                PAD = true;
        bool                SED = true;
        Real                tolerance = 1e-5;
        std::string         x_lower = "periodic";
        std::string         x_upper = "none";
        Real                xlim0 = 0.0;
        Real                xlim1 = 1.0;
        std::string         y_lower = "periodic";
        std::string         y_upper = "none";
        Real                ylim0 = 0.0;
        Real                ylim1 = 1.0;
        std::string         z_lower = "periodic";
        std::string         z_upper = "none";
        Real                zlim0 = 0.0;
        Real                zlim1 = 1.0;
    };

    SolverParams make_solver_params(const Config& config) {
        Integrator integrator;
        if      (config.integrator == "forward_euler") {integrator = Integrator::forward_euler;}
        else if (config.integrator == "rk4") {integrator = Integrator::rk4;}
        else if (config.integrator == "ssprk2") {integrator = Integrator::ssprk2;}
        else if (config.integrator == "ssprk3") {integrator = Integrator::ssprk3;}
        else {throw std::invalid_argument("unknown integrator: " + config.integrator);}

        Setup setup;
        if      (config.setup == "sinus") {setup = Setup::SINUS;}
        else if (config.setup == "square") {setup = Setup::SQUARE;}
        else {throw std::invalid_argument("setup must be sinus or square");}

        const auto parse_bc = [](const std::string& value) -> BC {
            if      (value == "none") {return BC::none;}
            else if (value == "periodic") {return BC::periodic;}
            else {throw std::invalid_argument("unknown boundary condition: " + value);}
        };

        const auto parse_fallback_cascade = [](const std::string& value) -> FallbackCascade {
            if      (value == "none") {return FallbackCascade::none;}
            else if (value == "full") {return FallbackCascade::full;}
            else if (value == "MUSCL") {return FallbackCascade::MUSCL;}
            else if (value == "MUSCL0") {return FallbackCascade::MUSCL0;}
            else {throw std::invalid_argument("unknown fallback cascade: " + value);}
        };

        const auto parse_MUSCL_limiter = [](const std::string& value) -> SlopeLimiter {
            if      (value == "none") {return SlopeLimiter::none;}
            else if (value == "minmod") {return SlopeLimiter::minmod;}
            else if (value == "moncen") {return SlopeLimiter::moncen;}
            else {throw std::invalid_argument("unknown MUSCL limiter: " + value);}
        };

        const HydroParams hydro{
            .gamma = config.gamma,
            .cfl = config.cfl
        };
        const UniformMeshParams mesh{
            .xlim0 = config.xlim0,
            .xlim1 = config.xlim1,
            .nx = static_cast<std::size_t>(config.nx),
            .ylim0 = config.ylim0,
            .ylim1 = config.ylim1,
            .ny = config.ny,
            .zlim0 = config.zlim0,
            .zlim1 = config.zlim1,
            .nz = config.nz
        };
        const BC_Params bc{
            .nghost = static_cast<std::size_t>(config.nghost),
            .x_lower = parse_bc(config.x_lower),
            .x_upper = parse_bc(config.x_upper),
            .y_lower = parse_bc(config.y_lower),
            .y_upper = parse_bc(config.y_upper),
            .z_lower = parse_bc(config.z_lower),
            .z_upper = parse_bc(config.z_upper)
        };
        const IOParams io{
            .fprint=config.fprint,
            .nout=config.nout,
            .output_dir=config.output_dir
        };
        const MOOD_Scheme scheme(
            config.p, 
            config.nrev, 
            parse_fallback_cascade(config.fallback_cascade), 
            parse_MUSCL_limiter(config.MUSCL_limiter),
            config.tolerance,
            config.NAD,
            config.PAD,
            config.SED
        );
        return {
            .hydro = hydro,
            .mesh = mesh,
            .bc = bc,
            .io = io,
            .scheme = scheme,
            .integrator = integrator,
            .setup = setup,
            .tstop = config.tstop,
            .max_steps = config.max_steps,
            .clip_last_dt = config.clip_last_dt
        };
    }

    std::string trim(const std::string& text) {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {return "";}
        return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
    }
}

SolverParams read_config(const std::string& path) {
    std::ifstream file(path);
    if (!file) {throw std::invalid_argument("cannot open config file: " + path);}

    Config config;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(file, line)) {
        ++line_number;
        line = trim(line.substr(0, line.find('#')));
        if (line.empty()) {continue;}
        const std::string location = path + ":" + std::to_string(line_number) + ": ";
        const auto equal = line.find('=');
        if (equal == std::string::npos) {throw std::invalid_argument(location + "expected key=value");}
        const std::string key = trim(line.substr(0, equal));
        std::istringstream input(line.substr(equal + 1));
        if      (key == "cfl") {input >> config.cfl;}
        else if (key == "clip_last_dt") {input >> std::boolalpha >> config.clip_last_dt;}
        else if (key == "fallback_cascade") {input >> config.fallback_cascade;}
        else if (key == "fprint") {input >> config.fprint;}
        else if (key == "gamma") {input >> config.gamma;}
        else if (key == "integrator") {input >> config.integrator;}
        else if (key == "max_steps") {input >> config.max_steps;}
        else if (key == "MUSCL_limiter") {input >> config.MUSCL_limiter;}
        else if (key == "NAD") {input >> std::boolalpha >> config.NAD;}
        else if (key == "nghost") {input >> config.nghost;}
        else if (key == "nout") {input >> config.nout;}
        else if (key == "nrev") {input >> config.nrev;}
        else if (key == "nx") {input >> config.nx;}
        else if (key == "ny") {input >> config.ny;}
        else if (key == "nz") {input >> config.nz;}
        else if (key == "output_dir") {input >> config.output_dir;}
        else if (key == "p") {input >> config.p;}
        else if (key == "PAD") {input >> std::boolalpha >> config.PAD;}
        else if (key == "SED") {input >> std::boolalpha >> config.SED;}
        else if (key == "setup") {input >> config.setup;}
        else if (key == "tolerance") {input >> config.tolerance;}
        else if (key == "tstop") {input >> config.tstop;}
        else if (key == "x_lower") {input >> config.x_lower;}
        else if (key == "x_upper") {input >> config.x_upper;}
        else if (key == "xlim0") {input >> config.xlim0;}
        else if (key == "xlim1") {input >> config.xlim1;}
        else if (key == "y_lower") {input >> config.y_lower;}
        else if (key == "y_upper") {input >> config.y_upper;}
        else if (key == "ylim0") {input >> config.ylim0;}
        else if (key == "ylim1") {input >> config.ylim1;}
        else if (key == "z_lower") {input >> config.z_lower;}
        else if (key == "z_upper") {input >> config.z_upper;}
        else if (key == "zlim0") {input >> config.zlim0;}
        else if (key == "zlim1") {input >> config.zlim1;}
        else {throw std::invalid_argument(location + "unknown config parameter: " + key);}
        char trailing;
        if (!input || (input >> trailing)) {throw std::invalid_argument(location + "invalid value for " + key);}
    }
    if (file.bad()) {throw std::invalid_argument("cannot read config file: " + path);}
    if (config.integrator.empty()) {throw std::invalid_argument(path + ": missing required key: integrator");}
    if (config.nghost == -1) {throw std::invalid_argument(path + ": missing required key: nghost");}
    if (config.nx == -1) {throw std::invalid_argument(path + ": missing required key: nx");}
    if (config.p == -1) {throw std::invalid_argument(path + ": missing required key: p");}
    if (config.setup.empty()) {throw std::invalid_argument(path + ": missing required key: setup");}
    if (config.tstop == -1.0) {throw std::invalid_argument(path + ": missing required key: tstop");}

    return make_solver_params(config);
}
