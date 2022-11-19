function hfit = dy_dt_eval(t,par,flist,y,M)

%evaluate dy/dt 
dy_dt = NaN*ones(length(t),3);
for p = 1:length(t)
    dy_dt(p,:) = ode_growth_linkingT3A(t(p),y(p,:),par,flist);
end
%linear fitting of dM/dt vs M --> the slope gives growth rate 
%becase dM/dt = lambda*M
hfit = fit(M(y),M(dy_dt),'poly1');

end