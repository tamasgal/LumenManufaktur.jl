"""
    Root(R, n, t)

Helper struct to find solution(s) to ``z`` of the square root expression:

```math
\begin{aligned}
ct(z=0) & = & z + n \\sqrt{z^2 + R^2}
\end{aligned}
```

where ``n = 1/\\cos(\\theta_{c})`` is the index of refraction.

# Arguments

- `R`: minimal distance of approach [m]
- `n`: index of refraction
- `t`: time at z = 0 [ns]

# Fields

- `most_upstream`: most upstream solution
- `most_downstream`: most downstream solution
- `isvalid`: validity of the solution
```
"""
struct Root
    most_upstream::Float64
    most_downstream::Float64
    isvalid::Bool

    function Root(R, n, t)
        a = n^2 - 1.0
        b = 2 * C * t
        c = R^2 * n^2 - C^2 * t^2

        q = b * b - 4 * a * c

        isvalid = false
        if q >= 0.0

            first = (-b - sqrt(q)) / (2a)
            second = (-b + sqrt(q)) / (2a)

            isvalid = C * t > second
        end
        new(first, second, isvalid)
    end
end
