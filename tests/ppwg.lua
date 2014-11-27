j=complex.I -- For convenience 

parameters = {}
for param = 0.1E12, 1.5E12, 0.05E12 do   parameters[#parameters+1] = param   end

-- Software Options
options = {    lambda_tgt = 1.4,
               output_dir = os.getenv("WORK").."/data_gr3d_base",
        update_lambda_tgt = false, 
         update_initspace = false,
           save_solutions = false 
           }

DIR = options['output_dir']
Edat = {DIR .. "/E0.dat", DIR .. "/E1.dat", DIR .. "/E2.dat"}
Ddat = {DIR .. "/D0.dat", DIR .. "/D1.dat", DIR .. "/D2.dat"}
Kdat = {DIR .. "/K0.dat", DIR .. "/K1.dat", DIR .. "/K2.dat"}

function p0(x)   return x^0   end
function p1(x)   return x^1   end
function p2(x)   return x^2   end

Efuncs = {p0, p1, p2}
Dfuncs = {p0, p1, p2}
Kfuncs = {p0, p1, p2}

