struct Distribution
	αa
	αb
	αc
	αd
	αe
	αf
	αg
	αh
end

function calculateVariable(point, BoxX, BoxY, X, Y)
	
	b = BoxY/2
	f = b

	d = BoxX/2
	h = d
	
	a = sqrt(h^2 + b^2)
	c = a
	e = a
	g = a

	B = point[1] - b
	F = Y - (point[1] + f)

	H = point[2] - h
	D = X - (point[2] + d)

	A = sqrt((H+h)^2 + (B+b)^2)
	C = sqrt((D+d)^2 + (B+b)^2)
	E = sqrt((D+d)^2 + (F+f)^2)
	G = sqrt((H+h)^2 + (F+f)^2)

	αa = A / a
	αb = B / b
	αc = C / c
	αd = D / d
	αe = E / e
	αf = F / f
	αg = G / g
	αh = H / h

	return Distribution(αa, αb, αc, αd, αe, αf, αg, αh)
end


function distributeArea(distribution::Distribution)

	
	
end
