function state_summary(state::ExtendedStefanState)
    return (
        t=state.t,
        max_speed=maximum(abs.(vec(state.speed_full))),
        n_frozen=count(state.frozen_mask),
    )
end
