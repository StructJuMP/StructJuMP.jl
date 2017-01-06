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

@variable(m, pmin[j] <= p[SCEN,j=NODE,TIME] <= pmax[j])
@variable(m, 0 <= dp[SCEN,filter(x->x.ltype=="a",LINK),TIME] <= 100)
@variable(m, 1 <= fin[SCEN,Link,TIME] <= 500)
@variable(m, 1 <= fout[SCEN,j=SUP,TIME] <= smax[j])
@variable(m, dem[SCEN,DEM,TIME])
@variable(m, 0 <= pow[SCEN,filter(x->x.ltype=="a",LINK),TIME] <= 3000)
@variable(m, 0 <= slack[SCEN,LINK,TIME,1:Nx])

@variable(m, 20 <= px[SCEN,LINK,TIME,DIS] <= 100)
@variable(m,  1 <= fx[SCEN,LINK,TIME,DIS] <= 500)

@variable(m, nu)
@variable(m, 0 <= phi[SCEN])
@variable(m, cvarcost)

for k in SCEN
    @constraint(m, cost[k] - nu <= phi[k])
end

@objective(m, (1-lambda)*mcost + lambda*cvarcost)

for k in SCEN, i in NODE, t in TIME
    @constraint(m, sum(fout[k,j,t] for j=LINK if j.lendoc==i)
                    + sum(s[k,j,t] for j=SUP if j.sloc==i)
                    - sum(fin[k,j,t] for j=LINK if j.lstartloc==i)
                    - sum(dem[k,j,t] for j=DEM if j.dloc==i) == 0)
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:(Nx-1)
    @constraint(m, (px[j,i,t+1,k]-px[j,i,t,k])/dt + c1[i]*(fx[j,i,t+1,k+1]-fx[j,i,t+1,k])/(dx[i]) == 0)
end

for j in SCEN, i in LINK, t in TIME
    @constraint(m, fx[j,i,t,1]  == fin[j,i,t])
    @constraint(m, fx[j,i,t,Nx] == fout[j,i,t])
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:(Nx-1)
    @constraint(m, (fx[j,i,t+1,k]-fx[j,i,t,k])/dt == -c2[i]*(px[j,i,t+1,k+1]-px[j,i,t+1,k])/dx[i] - slack[j,i,t+1,k])
end

for j in SCEN, i in LINK, t in TIMEm, k in 1:Nx
    @constraint(m, slack[j,i,t,k]*px[j,i,t,k] == c3[i]*fx[j,i,t,k]*fx[j,i,t,k])
end

for j in SCEN, i in filter(x->x.ltype=="p",LINK), t in TIME
    @constraint(m, px[j,i,t,1]  == p[j,lstartloc[i],t])
    @constraint(m, px[j,i,t,Nx] == p[j,lendloc[i],t])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in TIME
    @constraint(m, px[j,i,t,1]  == p[j,lstartloc[i],t] + dp[j,i,t])
    @constraint(m, px[j,i,t,Nx] == p[j,lendloc[i],t])
end

for k in SCEN, j in SUP, t in TIME
    @constraint(m, p[k,sloc[j],t] == pmin[sloc[j]])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in TIME
    @constraint(m, p[j,lstartloc[i],t] + dp[j,i,t] <= pmax[lstartloc[i]])
end

for j in SCEN, i in filter(x->x.ltype=="a",LINK), t in filter(x->(x<=TDEC),TIME)
    @constraint(m, dp[j,i,t] == dp[1,i,t])
end

for j in SCEN, i in LINK, t in filter(x->(x<=TDEC),TIME)
    @constraint(m, dem[j,i,t] == dem[1,i,t])
end

for j in SCEN, i in LINK, t in 0:0, k in 1:(Nx-1)
    @constraint(m, fx[j,i,t+1,k+1] - fx[j,i,t+1,k)] == 0)
    @constraint(m, 0 == -c2[i]*(px[j,i,t+1,k+1]-px[j,i,t+1,k])/dx[i] - slack[j,i,t+1,k])
end
