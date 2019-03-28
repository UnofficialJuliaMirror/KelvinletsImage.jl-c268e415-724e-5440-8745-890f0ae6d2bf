module variationApplier
    include("triangleInterpolator.jl")
    using ..KelvinletsObject, .triangleInterpolator, ..LinearAlgebra, ..Images, ..ProgressMeter
    export  __applyVariation__, applyAreaVariation__, __applyCloserAreaVariation__, __applyPerimeterBasedAreaVariation__, __applyEdgeBasedAreaVariation__, __applyPonderedAvgAreaVariation__

    function __applyVariation__(object::kelvinletsObject,
                                variationFunction::Function,
                                retardationFunction::Function,
                                heatmap::Bool
            )::Array{RGB{Float64}, 2}

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)

        for i=1:object.sizeY
            for j=1:object.sizeX
                Δ = variationFunction([i, j])

                dx1 = j - 1
                dx2 = object.sizeX - j
                dy1 = i - 1
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy) / object.sizeY
                x = 2(object.sizeX/2 - dx) / object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyAreaVariation__(object::kelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)

        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)

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

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyCloserAreaVariation__(object::kelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)

        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)
        for i=1:object.sizeY
            for j=1:object.sizeX

                df00 = norm([i, j] - [minY, minX])
                df01 = norm([i, j] - [minY, maxX])
                df10 = norm([i, j] - [maxY, minX])
                df11 = norm([i, j] - [maxY, maxX])

                closer = min(df00, df01, df10, df11)

                if closer == df00
                    Δ = variationFunction([i, j], [minY, minX])
                elseif closer == df01
                    Δ = variationFunction([i, j], [minY, maxX])
                elseif closer == df10
                    Δ = variationFunction([i, j], [maxY, minX])
                elseif closer == df11
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

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)

    end

    function __applyPerimeterBasedAreaVariation__(object::kelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                divided = [0., 0.]
                divisor = 0.

                for k=minY:maxY
                    for w=minX:maxX
                        d = norm([i, j] - [k, w]) + 0.0001
                        divided += (1/d) * [k, w]
                        divisor += (1/d)
                    end
                end

                newReference = divided/divisor
                Δ = variationFunction([i, j], round.(Int64, newReference))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyEdgeBasedAreaVariation__(object::kelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                d1 = norm([i, j] - [minY, minX]) + 1
                d2 = norm([i, j] - [minY, maxX]) + 1
                d3 = norm([i, j] - [maxY, minX]) + 1
                d4 = norm([i, j] - [maxY, maxX]) + 1

                newReference = ((1/d1) * [minY, minX] + (1/d2) * [minY, maxX] + (1/d3) * [maxY, minX] + (1/d4) * [maxY, maxX]) / ((1/d1) + (1/d2) + (1/d3) + (1/d4))
                Δ = variationFunction([i, j], round.(Int64, newReference))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])

                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)


                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyPonderedAvgAreaVariation__(object::kelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{Float64}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0.,0.,0.), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                d1 = norm([i, j] - [minY, minX]) + 1
                d2 = norm([i, j] - [minY, maxX]) + 1
                d3 = norm([i, j] - [maxY, minX]) + 1
                d4 = norm([i, j] - [maxY, maxX]) + 1

                deform1 = variationFunction([i,j], [minY, minX])
                deform2 = variationFunction([i,j], [minY, maxX])
                deform3 = variationFunction([i,j], [maxY, minX])
                deform4 = variationFunction([i,j], [maxY, maxX])

                Δ = ((1/d1) * deform1 + (1/d2) * deform1 + (1/d3) * deform1 + (1/d4) * deform1) /
                     ((1/d1) + (1/d2) + (1/d3) + (1/d4))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{Float64}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __interpolateVariation__(object::kelvinletsObject,
                                      allΔ::Array{Float64, 3}
            )::Array{RGB{Float64}, 2}

        interpImg = fill(RGB(0., 0., 0.), object.sizeY, object.sizeX)

        rasterize = function(A, B, C)
            colorA = object.image[A[1], A[2]]
            colorB = object.image[B[1], B[2]]
            colorC = object.image[C[1], C[2]]
            rasterizationBBOX(interpImg,
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

end
