function [i_set] = set_angle_calculator (alpha_h, alpha_f, i_w)

    alpha_w = i_w + alpha_f;
    eps_0 = 1;
    deps_dalfa = 0.3;
    eps = eps_0 + deps_dalfa * alpha_w;
    i_set = alpha_h - alpha_f + eps;

end