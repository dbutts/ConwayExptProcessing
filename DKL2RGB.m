function v_out = DKL2RGB(v_in, v_bkg, T_DKL2RGB, forward)

if nargin < 4
    forward = true;
end

if forward % convert DKL to RGB
v_out = v_bkg.* (T_DKL2RGB * v_in) + v_bkg;
else % convert RGB to DKL

    v_out = T_DKL2RGB \ ((v_in - v_bkg)./v_bkg)  ;
end

end
