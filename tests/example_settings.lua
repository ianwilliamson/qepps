j=complex.I --For convenience 

lambda_tgt = 55.5 - 6*j

parameters = {1E9, 1.2E9, 1.4E9, 1.6E9, 1.8E9, 2E9, 2.2E9} -- Explicit listing

--parameters = {} -- Building with a loop
--for param = 1E9, 2E9, 0.1E9 do   parameters[#parameters+1] = param   end


DIR = os.getenv("WORK") .. "/test_project"
Edat = {DIR .. "/E0.dat", DIR .. "/E1.dat", DIR .. "/E2.dat"}
Ddat = {DIR .. "/D0.dat", DIR .. "/D1.dat", DIR .. "/D2.dat"}
Kdat = {DIR .. "/K0.dat", DIR .. "/K1.dat", DIR .. "/K2.dat"}


local function p2(param)
   return param^2 
end

local function p1(param)
   return param^1 
end

local function p0(param)
   return param^0
end

Efuncs = {p0, p1, p2}
Dfuncs = {p0, p1, p2}
Kfuncs = {p0, p1, p2}

