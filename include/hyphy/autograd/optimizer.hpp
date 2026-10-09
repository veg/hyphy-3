#pragma once

#include "hyphy/autograd/var.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace hyphy::autograd {

class Adam {
public:
    std::vector<Var*> parameters;
    Scalar lr = 0.01;
    Scalar beta1 = 0.9;
    Scalar beta2 = 0.999;
    Scalar eps = 1e-8;

    std::vector<Scalar> m;
    std::vector<Scalar> v;
    size_t t = 0;

    Adam(std::vector<Var*> params, Scalar learning_rate = 0.01, Scalar b1 = 0.9, Scalar b2 = 0.999, Scalar epsilon = 1e-8)
        : parameters(std::move(params)), lr(learning_rate), beta1(b1), beta2(b2), eps(epsilon) {
        m.assign(parameters.size(), 0.0);
        v.assign(parameters.size(), 0.0);
    }

    void zero_grad() {
        for (auto* p : parameters) {
            if (p) p->zero_grad();
        }
    }

    void step() {
        t++;
        Scalar beta1_corr = 1.0 - std::pow(beta1, static_cast<Scalar>(t));
        Scalar beta2_corr = 1.0 - std::pow(beta2, static_cast<Scalar>(t));

        for (size_t i = 0; i < parameters.size(); ++i) {
            if (!parameters[i] || !parameters[i]->requires_grad()) continue;

            Scalar g = parameters[i]->grad();

            // Update biased first moment
            m[i] = beta1 * m[i] + (1.0 - beta1) * g;
            // Update biased second raw moment
            v[i] = beta2 * v[i] + (1.0 - beta2) * g * g;

            // Bias-corrected estimates
            Scalar m_hat = m[i] / beta1_corr;
            Scalar v_hat = v[i] / beta2_corr;

            // Gradient step
            Scalar new_val = parameters[i]->val() - lr * m_hat / (std::sqrt(v_hat) + eps);
            parameters[i]->set_val(new_val);
        }
    }
};

} // namespace hyphy::autograd
