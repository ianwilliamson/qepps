j=complex.I -- For convenience 

-- Parameter values
parameters = {}
for param = 4E12, 8E12, 0.25E12 do   parameters[#parameters+1] = param   end

-- Options
options = {}
options["lambda_tgt"] = 49.14-6.18*j --Initial target for eigenvalue
options["output_dir"] = os.getenv("WORK").."/data_gr3d_base" --Location to save solution vectors
options["output_log"] = options["output_dir"].."/output.txt" --File in which qepps (text) output will be saved
options["update_lambda_tgt"] = true --Update target eigenvalue from eigenvalue solved at previous parameter value
options["update_initspace"] = true --Update solver space from solution vector of previous parameter value
options["save_solutions"] = false --Save the solution vector for each parameter value
options["print_timing"] = true --At conclusion of parameter sweep, print timing

-- Scaling functions
function p0(x)   return x^0   end
function p1(x)   return x^1   end
function p2(x)   return x^2   end
-- Graphene surface conductivity calculation
function pS(x)
  gr_Ef = 0.4
  gr_gamma = 33.33
  pi = math.pi
  k0 = 2*pi*x/3E8
  omega_cm = 1e-2*k0/2/pi
  eSqDivh  = (1.602e-19)^2/(2*pi*1.05457148e-34)
  omegaF   = 1e4/1.24*gr_Ef
  sGR = complex.conj( 
        j*eSqDivh*(2*omegaF/(omega_cm+j*gr_gamma)+1/2*complex.log(complex.abs((omega_cm-2*omegaF)
        /(omega_cm+2*omegaF)))-j*pi/2*( 0 )) 
        )
  return x*sGR
end

-- Data files
matricies = {E={},D={},K={}}
matricies.E.data = {options["output_dir"].."/E2.dat"}
matricies.E.func = {p2}
matricies.D.data = {options["output_dir"].."/D1.dat"}
matricies.D.func = {p1}
matricies.K.data = {options["output_dir"].."/K0.dat",options["output_dir"].."/K2.dat",options["output_dir"].."/Ks.dat"}
matricies.K.func = {p0,p2,pS}
