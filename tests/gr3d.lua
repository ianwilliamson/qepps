j=complex.I -- For convenience 

-- Parameter values to sweep
parameters = {}
for param = 3E12, 15E12, 0.25E12 do   parameters[#parameters+1] = param   end

-- Software Options
options = {    lambda_tgt = 49.14-6.18*j, --Initial target for eigenvalue
               output_dir = os.getenv("WORK").."/data_gr3d_base", --
        update_lambda_tgt = true, --Should the eigenvalue target be updated after each parameter value
         update_initspace = false, --For each param value, should the solver's init vector space be updated  with the solution from the previous param value
           save_solutions = false --Should the solution vectors be saved for every parameter value
           }

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
matricies = {}
matricies["E"] = {options["output_dir"].."/E2.dat",p2}
matricies["D"] = {options["output_dir"].."/D1.dat",p1}
matricies["K"] = {options["output_dir"].."/K0.dat",p0,
                  options["output_dir"].."/K2.dat",p2,
                  options["output_dir"].."/Ks.dat",pS}

