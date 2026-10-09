#pragma once

#include "hyphy/core/types.hpp"
#include <string>
#include <vector>
#include <span>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <iostream>

namespace hyphy::core {

enum class ParamType {
    FREE,       // Standard free parameter in x
    FIXED,      // Fixed constant value (not in x)
    TIED,       // Shared with another parameter (points to index in x)
    LINEAR,     // theta = scale * x[k] + offset
    POWER       // theta = (x[k]) ^ power
};

struct ParamConstraint {
    ParamType type = ParamType::FREE;
    size_t opt_index = 0;       // Index into optimization vector x
    Scalar fixed_val = 1.0;     // Value if type == FIXED
    Scalar scale = 1.0;         // For LINEAR
    Scalar offset = 0.0;        // For LINEAR
    Scalar power = 1.0;         // For POWER
    std::string name;

    Scalar evaluate(std::span<const Scalar> x) const {
        switch (type) {
            case ParamType::FREE:
            case ParamType::TIED:
                return x[opt_index];
            case ParamType::FIXED:
                return fixed_val;
            case ParamType::LINEAR:
                return scale * x[opt_index] + offset;
            case ParamType::POWER:
                return std::pow(x[opt_index], power);
        }
        return fixed_val;
    }

    // Derivative d(theta) / d(x_k)
    Scalar derivative(std::span<const Scalar> x) const {
        switch (type) {
            case ParamType::FREE:
            case ParamType::TIED:
                return 1.0;
            case ParamType::FIXED:
                return 0.0;
            case ParamType::LINEAR:
                return scale;
            case ParamType::POWER:
                return power * std::pow(x[opt_index], power - 1.0);
        }
        return 0.0;
    }
};

class ParameterMap {
public:
    // Definition of free optimizer parameters x
    struct FreeParamInfo {
        std::string name;
        Scalar init_val = 1.0;
        Scalar lower_bound = 1e-6;
        Scalar upper_bound = 1e4;
    };

    std::vector<FreeParamInfo> free_params;
    std::vector<ParamConstraint> model_params;
    std::unordered_map<std::string, size_t> model_param_names;

    size_t num_free_params() const {
        return free_params.size();
    }

    size_t num_model_params() const {
        return model_params.size();
    }

    // Add a free optimizer parameter
    size_t add_free_param(const std::string& name, Scalar init_val, Scalar lb = 1e-6, Scalar ub = 1e4) {
        size_t idx = free_params.size();
        free_params.push_back({name, init_val, lb, ub});
        return idx;
    }

    // Register a model parameter
    size_t add_model_param(const std::string& name, const ParamConstraint& constraint) {
        size_t idx = model_params.size();
        ParamConstraint c = constraint;
        c.name = name;
        model_params.push_back(c);
        model_param_names[name] = idx;
        return idx;
    }

    // Map vector x to model parameters theta
    void evaluate_model_params(std::span<const Scalar> x, std::span<Scalar> theta_out) const {
        if (theta_out.size() != model_params.size()) {
            throw std::runtime_error("theta_out size does not match model_params count");
        }
        for (size_t i = 0; i < model_params.size(); ++i) {
            theta_out[i] = model_params[i].evaluate(x);
        }
    }

    // Backpropagate gradients: grad_x = sum_j (grad_theta_j * d_theta_j/d_x_k)
    void backprop_gradients(
        std::span<const Scalar> x,
        std::span<const Scalar> grad_theta,
        std::span<Scalar> grad_x_out
    ) const {
        std::fill(grad_x_out.begin(), grad_x_out.end(), 0.0);
        for (size_t j = 0; j < model_params.size(); ++j) {
            const auto& c = model_params[j];
            if (c.type == ParamType::FIXED) continue;
            Scalar d_theta = c.derivative(x);
            grad_x_out[c.opt_index] += grad_theta[j] * d_theta;
        }
    }

    // Convenience: Pin parameter to fixed value (for Null model tests)
    void fix_param(size_t model_idx, Scalar val) {
        if (model_idx >= model_params.size()) return;
        model_params[model_idx].type = ParamType::FIXED;
        model_params[model_idx].fixed_val = val;
    }

    void fix_param(const std::string& name, Scalar val) {
        auto it = model_param_names.find(name);
        if (it != model_param_names.end()) {
            fix_param(it->second, val);
        }
    }
};

} // namespace hyphy::core
