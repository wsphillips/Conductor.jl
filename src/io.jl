# Return ODESystem pretty printing for our wrapper types
#Base.show(io::IO, ::MIME"text/plain", x::IonChannel) = Base.display(isbuilt(x) ? x.sys : x)
