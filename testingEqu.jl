include("src/KelvinletsImage.jl")
using ImageView, TestImages, Images, LinearAlgebra, Plots, .KelvinletsImage

img = convert(Array{RGB{Float64},2}, testimage("mandrill"))
grid = convert(Array{RGB{Float64},2}, load("testImages/grid.jpg"))
gridObj = kelvinletsObject(grid, 0.4, 1.)
obj = kelvinletsObject(img, 0.4, 1.)

new = grabRectangle(obj, [175 175; 250 250], [100., 100.], 50., false)

plot(plot(new), size=(512,512))
png("/home/ghaetinger/Desktop/allPixelsRef.png")
points = [100 100; 75 125; 100 150; 150 150; 175 125; 150 100]

new = grabRectangleALL(obj, [175 175; 200 200], [200., 200.], 250., false)
heat = grabRectangleALL(obj, [175 175; 200 200], [200., 200.], 250., true)
newGrid = grabRectangleALL(gridObj, [175 175; 200 200], [200., 200.], 250., false)
unrestr = load("/home/ghaetinger/Desktop/hexaUnrestr.png")
plot(plot(new), plot(heat), plot(newGrid), size=(1000, 1000))
png("/home/ghaetinger/Desktop/hexaPoly.png")
png("/home/gghaetinger/Desktop/200,200,250.png")

save("/home/ghaetinger/Desktop/hexaUnrestr.png", unrestr)
new = load("/home/ghaetinger/Desktop/hexaNew.png")

vid = KelvinletsImage.makeVideo(obj, KelvinletsImage.grabRectangle__new, [175 175; 250 250], [100., 100.], 50., 20)
imshow(vid)


points[size(points[1]), :]
size(points)[1]

[100 100; 200 200]

norm([100, 100] - [70, 125])
norm([70, 125] - [22, 165])
norm([22, 165] - [100, 100])


include("src/geometryAnalisys.jl")

a = [100, 100]
b = [100, 200]
c = [200, 200]
d = [200, 100]

plot([geometryAnalisys.reference(a, b, d, c, [150, i]) for i = 1:500])


retardationFunc = α -> (cos(π * α) + 1) / 2
retardationFunc = α -> (2/sqrt(2))*(cos(π * α) + 1)^(1/2) / 2
retardationFunc = α -> (1/(2^(1/4)))*(cos(π * α) + 1)^(1/4)
ls = []
for j = 1:500
    dx1 = j - 1
    dx2 = 500 - j
    dx = min(dx1, dx2)
    x = 2(500/2 - dx) / 500
    append!(ls, x)
end
plot(plot(1:500, retardationFunc.(ls), label="variation on the X axis"), plot(new), size=(1000, 500))
png("/home/ghaetinger/Desktop/newRetResult.png")

minx = 200
maxx = 300
mean = (200 + 300)/2
imgmax = 500

ls = []
for i = 0:imgmax
    if i < minx
        x = i/minx
        a = retardationFunc(x)
    elseif i > maxx
        x = (imgmax - i)/(imgmax - maxx)
        a = retardationFunc(x)
    else
        a = 0
    end
    append!(ls, a)
end

plot(ls, label="distance coeficient")
png("/home/ghaetinger/Desktop/distanceCoef.png")

x = 0:0.1:10

plot(x, retardationFunc.(x))

retardationFunc(1)

@show a = 10:-1:1
