module KelvinletsImage

    include("triangleInterpolator.jl")
    using Images, ProgressMeter, ImageView, LinearAlgebra
    export KelvinletsObject, grab, scale, pinch, grabRectangle, makeVideo
    
    struct KelvinletsObject
        a::Float64
        b::Float64
        c::Float64
        sizeX::Int64
        sizeY::Int64
        image::AbstractArray{RGB{N0f8}, 2}
        function KelvinletsObject(image::AbstractArray{RGB{N0f8}, 2},
                                  ν::Float64,
                                  μ::Float64
                )::KelvinletsObject

            a = 1 / (4pi * μ)
            b = a / (4(1 - ν))
            c = 2 / (3a- 2b)

            sizeY, sizeX = size(image)

            new(a, b, c, sizeX, sizeY, image)
        end
    end

    function __applyVariation__(object::KelvinletsObject,
                                pressurePoint::Array{Int64},
                                variationFunction::Function,
                                retardationFunction::Function
            )::Array{RGB{N0f8}, 2}

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        for i=1:object.sizeY
            for j=1:object.sizeX
                Δ = variationFunction([i, j])

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy)/object.sizeY
                x = 2(object.sizeX/2 - dx)/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function
            )::Array{RGB{N0f8}, 2}
        
        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        for i=1:object.sizeY
            for j=1:object.sizeX

                if i < minY && j < minX #Top-LEFT corner
                    Δ = variationFunction([i, j], [minY, minX])
                elseif i < minY && j >= minX && j <= maxX # Top
                    Δ = variationFunction([i, j], [minY, j])
                elseif i < minY && j > maxX #Top-RIGHT corner
                    Δ = variationFunction([i, j], [minY, maxX])
                elseif j < minX && i >= minY && i <= maxY #LEFT
                    Δ = variationFunction([i, j], [i, minX])
                elseif j >= minX && j <= maxX && i >= minY && i <= maxY # MIDDLE
                    Δ = variationFunction([i, j], [i, j])
                elseif j > maxX && i >= minY && i <= maxY #RIGHT
                    Δ = variationFunction([i, j], [i, maxX])
                elseif i > maxY && j < minX #Lower-LEFT corner
                    Δ = variationFunction([i, j], [maxY, minX])
                elseif i > maxY && j >= minX && j <= maxX #BASE
                    Δ = variationFunction([i, j], [maxY, j])
                elseif i >  maxY && j > maxX #Lower-RIGHT corner
                    Δ = variationFunction([i, j], [maxY, maxX])
                end

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy)/object.sizeY
                x = 2(object.sizeX/2 - dx)/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __interpolateVariation__(object::KelvinletsObject,
                                      allΔ::Array{Float64, 3}
            )::Array{RGB{N0f8}, 2}

        interpImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)
        
        rasterize = function(A, B, C)
            colorA = object.image[A[1], A[2]]
            colorB = object.image[B[1], B[2]]
            colorC = object.image[C[1], C[2]]
            triangleInterpolator.rasterizationBBOX(interpImg,
                            [allΔ[A[1], A[2], 1], allΔ[A[1], A[2], 2]],
                            [allΔ[B[1], B[2], 1], allΔ[B[1], B[2], 2]],
                            [allΔ[C[1], C[2], 1], allΔ[C[1], C[2], 2]],
                            colorA, colorB, colorC
            )
        end
        
        for i=1:object.sizeY-1
          for j=1:object.sizeX-1
            rasterize([i, j], [i, j+1], [i+1, j+1])
            rasterize([i, j], [i+1, j], [i+1, j+1])
          end
        end
        return interpImg
    end

    """
        grab(object::KelvinletsObject, x0::Array{Int64}, force::Array{Float64}, ϵ::Float64)

    grabs a point in an image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a force vector *force*
    and a brush size *ϵ*
    # Example:
    ```julia-repl
        grab(obj, [100, 100], [100., 0.], 50.)
    ```
    """    
    function grab(object::KelvinletsObject,
                  x0::Array{Int64},
                  force::Array{Float64},
                  ϵ::Float64
            )::Array{RGB{N0f8}, 2}

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
        return __applyVariation__(object, x0, grabFunc, retardationFunc)
    end

    """
        scale(object::KelvinletsObject, x0::Array{Int64}, scale::Float64, ϵ::Float64)

    Scales a point on an image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a scale alpha *scale*
    and a brush size *ϵ*
    # Example:
    ```julia-repl
        julia> scale(obj, [100, 100], 50., 50.)
    ```
    """
    function scale(object::KelvinletsObject,
                   x0::Array{Int64},
                   force::Float64,
                   ϵ::Float64
            )::Array{RGB{N0f8}, 2}

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
        return __applyVariation__(object, x0, scaleFunc, retardationFunc)
    end
    
    """
         pinch(object::KelvinletsObject, x0::Array{Int64}, force::Array{Float64, 2}, ϵ::Float64)

    Pinches the image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a force matrix *force*
    and a brush size *ϵ*
    # Example:
    ```julia-repl
        julia> pinch(obj, [100, 100], [10000. 0. ; 10000. 0.], 50.)
    ```
    """
    function pinch(object::KelvinletsObject,
                   x0::Array{Int64},
                   force::Array{Float64, 2},
                   ϵ::Float64
            )::Array{RGB{N0f8}, 2}

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
        return __applyVariation__(object, x0, pinchFunc, retardationFunc)
    end

    """
        grabRectangle(object::KelvinletsObject, points::Array{Int64, 2}, force::Array{Float64}, ϵ::Float64)

    Grabs the image using a square brush given a KelvinletsObject *obj*,
    a given matrix representation of a rectangle *points* in the following format --> [minY, minX ; maxY, maxX],
    a force vector *force*,
    and a given brush size *ϵ* (to calculate te variation of outside pixels)
    # Example:
    ```julia-repl
        julia> grabRectangle(obj, [100 100 ; 250 250], [100., 0.], 50.)
    ```
    """
    function grabRectangle(object::KelvinletsObject,
                           points::Array{Int64, 2},
                           force::Array{Float64},
                           ϵ::Float64
            )::Array{RGB{N0f8}, 2}

        grabFunc = function(x::Array{Int64},
                            x0::Array{Int64}
                    )
        
            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyAreaVariation__(object, points, grabFunc, retardationFunc)

    end

    """
        makeVideo(object::KelvinletsObject, kelvinletsFunction::Function, x0::Array{Int64}, force, ϵ::Float64, frames::Int64)

    Computes a video for a given KelvinletsImage deformation function *KelvinletsFunction*, 
    using a *KelvinletsObject* as reference,
    a *x0* point of force application,
    a *force* matrix/vector (depending on the function),
    a *ϵ* brush size
    and a number of frames *frames*
    # Example
    ```julia-repl
        makeVideo(object, grab, [100, 100], [200., 0.], 70., 20)
    ```
    """
    function makeVideo(object::KelvinletsObject,
                       kelvinletsFunction::Function,
                       x0::Array{Int64},
                       force,
                       ϵ::Float64,
                       frames::Int64
        )::Array{RGB{N0f8},3}
        
        if typeof(force) == Float64
            var = range(0, stop=force, length=frames)
        else    
            var = range(fill(0, size(force)), stop=force, length=frames)    
        end
        
        video = Array{RGB{N0f8}}(undef, object.sizeY, object.sizeX, frames)
        @showprogress for i=1:frames
            video[:,:,i] = kelvinletsFunction(object, x0, var[i], ϵ)
        end
        imshow(video)
        return video
    end
end
