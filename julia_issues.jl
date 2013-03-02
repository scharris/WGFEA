
# Ranges of short types are apparantly promoted up to larger types, bug?  Here Uint8 becomes Uint32 so no method matches.
uint8Fun(x::Uint8) = 1
l = [uint8Fun(k, dim) for k=zero(Uint8):uint8(10)]


# Would be very convenient if let variables could depend on previous ones in the same list of definitions:
let x = 1, y = 2; x + y end

# Can function taking two or more tuples be defined with destructuring?  The following gives an error.  (Otherwise, how to access tuple parts?)
distance((x1,y1),(x2,y2)) = sqrt((x2-x1)^2 + (y2-y1)^2)

