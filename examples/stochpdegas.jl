using StructJuMP

TIME = 1:Nt
TIMEm = 1:(Nt-1)
DIS = 1:Nx
SCEN = 1:S

epsil = 0.025
z = 0.8
rhon = 0.72
R = 8314.
M = 18.
Tgas = 293.15
Cp = 2.34
Cv = 1.85
gam = Cp / Cv
om = (gam-1) / gam

ce = 0.1
cd = 1e6
cT = 1e6
cs = 0.

dt = Tf / Nt

m = StructuredModel()

@defVar(m, pmin[j] <= p[SCEN,j=NODE,TIME] <= pmax[j])
@defVar(m, 0 <= dp[SCEN,filter(x->x.ltype=="a",LINK),TIME] <= 100)
@defVar(m, 1 <= fin[SCEN,Link,TIME] <= 500)
@defVar(m, 1 <= fout[SCEN,j=SUP,TIME] <= smax[j])
@defVar(m, dem[SCEN,DEM,TIME])
@defVar(m, 0 <= pow[SCEN,filter(x->x.ltype=="a",LINK),TIME] <= 3000)
@defVar(m, 0 <= slack[SCEN,LINK,TIME,1:Nx])

@defVar(m, 20 <= px[SCEN,LINK,TIME,DIS] <= 100)
@defVar(m,  1 <= fx[SCEN,LINK,TIME,DIS] <= 500)

@defVar(m, nu)
@defVar(m, 0 <= phi[SCEN])
@defVar(m, cvarcost)

for k in SCEN
    @addConstraint(m, cost[k] - nu <= phi[k])
end

@setObjective(m, (1-lambda)*mcost + lambda*cvarcost)

for k in SCEN, i in NODE, t in TIME
    @addConstraint(m, sum{fout[k,j,t], j=LINK; j.lendoc==i} 
                    + sum{s[k,j,t], j=SUP; j.sloc==i}
                    - sum{fin[k,j,t], j=LINK; j.lstartloc==i}
                    - sum{dem[k,j,t], j=DEM; j.dloc==i} == 0)
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:(Nx-1)
    @addConstraint(m, (px[j,i,t+1,k]-px[j,i,t,k])/dt + c1[i]*(fx[j,i,t+1,k+1]-fx[j,i,t+1,k])/(dx[i]) == 0)
end

for j in SCEN, i in LINK, t in TIME
    @addConstraint(m, fx[j,i,t,1]  == fin[j,i,t])
    @addConstraint(m, fx[j,i,t,Nx] == fout[j,i,t])
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:(Nx-1)
    @addConstraint(m, (fx[j,i,t+1,k]-fx[j,i,t,k])/dt == -c2[i]*(px[j,i,t+1,k+1]-px[j,i,t+1,k])/dx[i] - slack[j,i,t+1,k])
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:Nx
    @addConstraint(m, slack[j,i,t,k]*px[j,i,t,k] == c3[i]*fx[j,i,t,k]*fx[j,i,t,k])
end

for j in SCEN, i in filter(x->x.ltype=="p",LINK), t in TIME
    @addConstraint(m, px[j,i,t,1]  == p[j,lstartloc[i],t])
    @addConstraint(m, px[j,i,t,Nx] == p[j,lendloc[i],t])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in TIME
    @addConstraint(m, px[j,i,t,1]  == p[j,lstartloc[i],t] + dp[j,i,t])
    @addConstraint(m, px[j,i,t,Nx] == p[j,lendloc[i],t])
end

for k in SCEN, j in SUP, t in TIME
    @addConstraint(m, p[k,sloc[j],t] == pmin[sloc[j]])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in TIME
    @addConstraint(m, p[j,lstartloc[i],t] + dp[j,i,t] <= pmax[lstartloc[i]])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in filter(x->(x<=TDEC),TIME)
    @addConstraint(m, dp[j,i,t] == dp[1,i,t])
end

for j in SCEN, i in LINK, t in filter(x->(x<=TDEC),TIME)
    @addConstraint(m, dem[j,i,t] == dem[1,i,t])
end

for j in SCEN, i in LINK, t in 0:0, k in 1:(Nx-1)
    @addConstraint(m, fx[j,i,t+1,k+1] - fx[j,i,t+1,k)] == 0)
    @addConstraint(m, 0 == -c2[i]*(px[j,i,t+1,k+1]-px[j,i,t+1,k])/dx[i] - slack[j,i,t+1,k])
end
