x = range(0,2,15)

y = sin.(x)

z = Δ * y

w = similar(x)


dx = 1/14
for i in 2:n-1
    w[i]  = Lap2(y,i)/dx^2
end

w[1] = (y[2]-y[1])/dx^2
w[end] = (y[end-1]-y[end])/dx^2