module KelvinletsObject
    using ..Images
    export kelvinletsObject
    """
        kelvinletsObject(image::AbstractArray{RGB{Float64}, 2}, ν::Float64, μ::Float64)

    Initializes KelvinletsObject for a given image *image*,
    poisson ratio *ν*
    and elastic shear modulus *μ*

    # Example:
    ```julia-repl
    julia> object = kelvinletsObject(image, 0.4, 1.)
    ```
    """
    mutable struct kelvinletsObject
        a::Float64
        b::Float64
        c::Float64
        sizeX::Int64
        sizeY::Int64
        image::AbstractArray{RGB{Float64}, 2}
        function kelvinletsObject(image::AbstractArray{RGB{Float64}, 2},
                                  ν::Float64,
                                  μ::Float64
                )::kelvinletsObject

            a = 1 / (4pi * μ)
            b = a / (4(1 - ν))
            c = 2 / (3a- 2b)

            sizeY, sizeX = size(image)

            new(a, b, c, sizeX, sizeY, image)
        end
    end
end
