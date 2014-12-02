j=complex.I -- For convenience 

parameters = {}
for param = 0.1E12, 1.5E12, 0.05E12 do   parameters[#parameters+1] = param   end

-- Software Options
options = {    lambda_tgt = 1.4,
               output_dir = os.getenv("WORK").."/data_ppwg_base",
        update_lambda_tgt = false, 
         update_initspace = false,
           save_solutions = false,
             print_timing = true
           }

function p0(x)   return x^0   end
function p1(x)   return x^1   end
function p2(x)   return x^2   end

-- Data files
matricies = {}
matricies["E"] = {options["output_dir"].."/E2.dat",p2}
matricies["D"] = {options["output_dir"].."/D1.dat",p1}
matricies["K"] = {options["output_dir"].."/K0.dat",p0,
                  options["output_dir"].."/K2.dat",p2}
