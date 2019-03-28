module geometryAnalisys
    using ..LinearAlgebra
    export reference, polygonReference

    function triangleArea(a, b, x)
        ab = norm(a - b)
        bx = norm(b - x)
        xa = norm(x - a)
        s = (ab+bx+xa)/2
        toSq = s*(s-ab)*(s-bx)*(s-xa)
        if toSq < 0
            return 0.
        end
        return sqrt(toSq)
    end

    function reference(a, b, d, c, x, X, Y)
        #
        ax = a - x
        bx = b - x
        cx = c - x
        dx = d - x

        mainDiag = a - c

        Aab = triangleArea(a, b, x)
        Abc = triangleArea(b, c, x)
        Acd = triangleArea(c, d, x)
        Ada = triangleArea(d, a, x)

        At = abs(mainDiag[1]) * abs(mainDiag[2])

        da = norm(ax) + 0.00001
        db = norm(bx) + 0.00001
        dc = norm(cx) + 0.00001
        dd = norm(dx) + 0.00001

        # if(x[1] >= a[1] && x[1] <= c[1] && x[1] <= c[2] && x[2] >= a[2])
        #     if At !=
        # end

        return ((At/(Aab + Abc + Acd + Ada))-1)^2 *
                   (x -
                    ((1/da) * a + (1/db) * b + (1/dc) * c + (1/dd) * d) /
                    ((1/da) + (1/db) + (1/dc) + (1/dd))
                   )
        # minx = a[2]
        # maxx = c[2]
        # miny = a[1]
        # maxy = c[1]
        # i = x[1]
        # j = x[2]
        #
        # meanx = (maxx + minx)/2
        # meany = (maxy + miny)/2
        #
        # retardationFunc = α -> (cos(π * α) + 1) / 2
        #
        # if j < minx
        #     aX = j/minx
        #     aX = retardationFunc(aX)
        # elseif j > maxx
        #     aX = (X - j)/(X - maxx)
        #     aX = retardationFunc(aX)
        # else
        #     aX = 0
        # end
        #
        # if i < miny
        #     aY = i/maxy
        #     aY = retardationFunc(aY)
        # elseif i > maxx
        #     aY = (Y - j)/(Y - maxy)
        #     aY = retardationFunc(aY)
        # else
        #     aY = 0
        # end
        #
        #
        # return  (aX + aY)/2 * (x -
        #             ((1/da) * a + (1/db) * b + (1/dc) * c + (1/dd) * d) /
        #             ((1/da) + (1/db) + (1/dc) + (1/dd))
        #            )

    end

    function polygonReference(pointList, x)
        totalArea = 0.
        for i=2:size(pointList)[1] - 1
            totalArea += triangleArea(pointList[i, :], pointList[i+1, :], pointList[1, :])
        end

        sumOfAreas = 0.
        for i=1:size(pointList)[1] - 1
            sumOfAreas += triangleArea(pointList[i, :], pointList[i+1, :], x)
        end
        sumOfAreas += triangleArea(pointList[size(pointList)[1], :], pointList[1, :], x)

        divided = [0., 0.]
        divisor = 0.
        for i=1:size(pointList)[1]
            d = norm(pointList[i, :] - x) + 0.00000001
            divided += (1/d) * pointList[i, :]
            divisor += (1/d)
        end

        return (log(10, totalArea/sumOfAreas))^2 *
                   (x - (divided/divisor))
    end

end
