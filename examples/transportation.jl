using StructJuMP

m = StructuredModel()

factories = 1:3
centers = 1:5
scenarios = ["lo", "mid", "hi"]

cap = [500, 450, 650]

demand = [150 160 170;
          100 120 135;
          250 270 300;
          300 325 350;
          600 700 800]

prob = [0.25, 0.5, 0.25]

transcost = [2.49 5.21 3.76 4.85 2.07;
             1.46 2.54 1.83 1.86 4.76;
             3.26 3.08 2.60 3.76 4.45]

prodcost = 14
price = 24
wastecost = 4

@variable(m, ship[factories,centers] >= 0)
@variable(m, product[factories] >= 0)
@variable(m, sales[centers] >= 0)
@variable(m, waste[centers] >= 0)
@variable(m, profit)

@variable(m, received[centers] >= 0)
@objective(m, Max, -sum(transcost[i,j]*ship[i,j] for i=factories, j=centers) + sum(prodcost*product[i] for i=factories))
for j in centers
    @constraint(m, received[j] == sum(ship[i,j] for i=factories))
end
for i in factories
    @constraint(m, product[i] == sum(ship[i,j] for j=centers))
end

for (s, elem) in enumerate(scenarios)
    bl = StructuredModel(parent=m, id=s)
    @variable(bl, 0 <= salesw[i=centers] <= demand[i,s])
    @variable(bl, wastew[centers] >= 0)
    @objective(bl, :Max, sum(price*prob[s]*salesw[j] for j=centers) - sum(wastecost*prob[s]*wastew[j] for j=centers))
    for j in centers
        @constraint(bl, received[j] == salesw[j]+wastew[j])
    end
end
