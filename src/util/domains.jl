
iscompact(d::Domain) = false
iscompact(d::AbstractInterval) = true
iscompact(d::PeriodicInterval) = true
