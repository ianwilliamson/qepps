j=complex.I -- For convenience 

-- Parameter values
parameters = {}
for param = 0.1E12, 1.5E12, 0.05E12 do   parameters[#parameters+1] = param   end

-- Options
options = {}
options["lambda_tgt"] = 1.4 --Initial target for eigenvalue
options["output_dir"] = os.getenv("WORK").."/data_ppwg_base" --Location to save solution vectors
options["output_log"] = options["output_dir"].."/output.txt" --File in which qepps (text) output will be saved
options["update_lambda_tgt"] = false --Update target eigenvalue from eigenvalue solved at previous parameter value
options["update_initspace"] = false --Update solver space from solution vector of previous parameter value
options["save_solutions"] = false --Save the solution vector for each parameter value
options["print_timing"] = true --At conclusion of parameter sweep, print timing

-- Scaling functions
function p0(x)   return x^0   end
function p1(x)   return x^1   end
function p2(x)   return x^2   end

-- Data files
matricies = {E={},D={},K={}}
matricies.E.data = {options["output_dir"].."/E2.dat"}
matricies.E.func = {p2}
matricies.D.data = {options["output_dir"].."/D1.dat"}
matricies.D.func = {p1}
matricies.K.data = {options["output_dir"].."/K0.dat",options["output_dir"].."/K2.dat"}
matricies.K.func = {p0,p2}
