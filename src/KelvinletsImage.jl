module KelvinletsImage
    using Images, ProgressMeter, ImageView, LinearAlgebra
    include("geometryAnalisys.jl"); include("KelvinletsObject.jl"); include("variationApplier.jl")
    using .KelvinletsObject, .variationApplier, .geometryAnalisys
    export kelvinletsObject, grab, scale, pinch, grabRectangle, grabPolygon, makeVideo

    """
        grab(object::kelvinletsObject, x0::Array{Int64}, force::Array{Float64}, ϵ::Float64, heatmap::Bool)

        grabs a point in an image given a kelvinletsObject *obj*,
        a pressure point *x0*,
        a force vector *force*,
        a brush size *ϵ*
        and whether you would like to generate the heatmap of the deformation
        # Example:
        ```julia-repl
            julia> newImage = grab(obj, [100, 100], [100., 0.], 50.)
        ```
    """
    function grab(object::kelvinletsObject,
                  x0::Array{Int64},
                  force::Array{Float64},
                  ϵ::Float64,
                  heatmap::Bool
            )::Array{RGB{Float64}, 2}

        grabFunc = function(x::Array{Int64})

            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, grabFunc, retardationFunc, heatmap)
    end

    """
        scale(object::kelvinletsObject, x0::Array{Int64}, scale::Float64, ϵ::Float64, heatmap::Bool)

        Scales a point on an image given a kelvinletsObject *obj*,
        a pressure point *x0*,
        a scale alpha *scale*,
        a brush size *ϵ*,
        and whether you would like to generate the heatmap of the deformation
        # Example:
        ```julia-repl
            julia> newImage = scale(obj, [100, 100], 50., 50.)
        ```
    """
    function scale(object::kelvinletsObject,
                   x0::Array{Int64},
                   force::Float64,
                   ϵ::Float64,
                   heatmap::Bool
            )::Array{RGB{Float64}, 2}

        scaleFunc = function(x::Array{Int64})

            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)

            return (2 * object.b - object.a) *
                   ( (1 / rϵ^2) +
                   ((ϵ^2)) / (2 * (rϵ^4))) *
                   (force * r)
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, scaleFunc, retardationFunc, heatmap)
    end

    """
         pinch(object::kelvinletsObject, x0::Array{Int64}, force::Array{Float64, 2}, ϵ::Float64, heatmap::Bool)

        Pinches the image given a kelvinletsObject *obj*,
        a pressure point *x0*,
        a force matrix *force*,
        a brush size *ϵ*,
        and whether you would like to generate the heatmap of the deformation
        # Example:
        ```julia-repl
            julia> newImage = pinch(obj, [100, 100], [10000. 0. ; 10000. 0.], 50.)
        ```
    """
    function pinch(object::kelvinletsObject,
                   x0::Array{Int64},
                   force::Array{Float64, 2},
                   ϵ::Float64,
                   heatmap::Bool
            )::Array{RGB{Float64}, 2}

        pinchFunc = function(x::Array{Int64})

            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            return  -2 * object.a * ((1 / rϵ^2) +
                    (ϵ^2 / rϵ^4)) * force * r +
                    4 * object.b * ((1 / rϵ^2) * force -
                    (1 / rϵ^4) * (r' * force * r) * I) * r
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, pinchFunc, retardationFunc, heatmap)
    end

    """
        grabRectangle(object::kelvinletsObject, points::Array{Int64, 2}, force::Array{Float64}, ϵ::Float64, heatmap::Bool)

        Grabs the image using a square brush given a kelvinletsObject *obj*,
        a given matrix representation of a rectangle *points* in the following format --> [minY, minX ; maxY, maxX],
        a force vector *force*,
        a given brush size *ϵ* (to calculate te variation of outside pixels),
        and whether you would like to generate the heatmap of the deformation
        # Example:
        ```julia-repl
            julia> newImage = grabRectangle(obj, [100 100 ; 250 250], [100., 0.], 50.)
        ```
    """
    function grabRectangle(object::kelvinletsObject,
                           points::Array{Int64, 2},
                           force::Array{Float64},
                           ϵ::Float64,
                           heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        a = [minY, minX]
        b = [minY, maxX]
        c = [maxY, minX]
        d = [maxY, maxX]

        grabFunc = function(x::Array{Int64})

            r = reference(a, b, c, d, x, object.sizeX, object.sizeY)
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        # retardationFunc = α -> (2/sqrt(2))*(cos(π * α) + 1)^(1/2) / 2
        # retardationFunc = α -> (1/(2^(1/4)))*(cos(π * α) + 1)^(1/4)
        return __applyVariation__(object, grabFunc, retardationFunc, heatmap)
    end

    """
        grabPolygon(object::kelvinletsObject, points::Array{Int64, 2}, force::Array{Float64}, ϵ::Float64, heatmap::Bool)

        Grabs the image using a polygon brush given a kelvinletsObject *obj*,
        a given matrix representation of a polygon *points* in the following format --> [y1 x1 ; y2 x2; y3 x3; ...],
        a force vector *force*,
        a given brush size *ϵ* (to calculate te variation of outside pixels),
        and whether you would like to generate the heatmap of the deformation
        # Example:
        ```julia-repl
            julia> newImage = grabPolygon(obj, [100 100 ; 250 250; 100 250], [100., 0.], 50.)
        ```
    """
    function grabPolygon(object::kelvinletsObject,
                         points::Array{Int64, 2},
                         force::Array{Float64},
                         ϵ::Float64,
                         heatmap::Bool
            )::Array{RGB{Float64}, 2}

        grabFunc = function(x::Array{Int64})

            r = polygonReference(points, x)
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        # retardationFunc = α -> (2/sqrt(2))*(cos(π * α) + 1)^(1/2) / 2
        # retardationFunc = α -> (1/(2^(1/4)))*(cos(π * α) + 1)^(1/4)
        return __applyVariation__(object, grabFunc, retardationFunc, heatmap)
    end

    """
        makeVideo(object::kelvinletsObject, kelvinletsFunction::Function, x0::Array{Int64}, force, ϵ::Float64, frames::Int64)

        Computes a video for a given KelvinletsImage deformation function *KelvinletsFunction*,
        using a *kelvinletsObject* as reference,
        a *x0* point of force application,
        a *force* matrix/vector (depending on the function),
        a *ϵ* brush size
        and a number of frames *frames*
        # Example
        ```julia-repl
            julia> video = (object, grab, [100, 100], [200., 0.], 70., 20)
        ```
    """
    function makeVideo(object::kelvinletsObject,
                       kelvinletsFunction::Function,
                       x0::Array{Int64},
                       force,
                       ϵ::Float64,
                       frames::Int64
        )::Array{RGB{Float64},3}

        if typeof(force) == Float64
            var = range(0, stop=force, length=frames)
        else
            var = range(fill(0, size(force)), stop=force, length=frames)
        end

        video = Array{RGB{Float64}}(undef, object.sizeY, object.sizeX, frames)
        @showprogress for i=1:frames
            video[:,:,i] = kelvinletsFunction(object, x0, var[i], ϵ, false)
        end
        imshow(video)
        return video
    end

end
