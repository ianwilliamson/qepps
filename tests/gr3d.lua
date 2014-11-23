j=complex.I --For convenience 
lambda_tgt = 44-6*j
parameters = {1, 2, 3} -- Explicit listing

DIR = os.getenv("HOME") .. "/projects/qepps/tests/comsol_gr3d"
Edat = {DIR .. "/E.dat"}
Ddat = {DIR .. "/D.dat"}
Kdat = {DIR .. "/K.dat"}

local function p2(param)
   return param^2 
end

local function p1(param)
   return param^1 
end

local function p0(param)
   return param^0
end

Efuncs = {p1}
Dfuncs = {p1}
Kfuncs = {p1}

