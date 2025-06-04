using PlotlyJS
include("ControlUN.jl")
G = tf([1],[1,2,0])
plot(rlocus(G))
using PlotlyJS

using PlotlyJS

using PlotlyJS

traces = GenericTrace[]
for i in 1:4
    trace = scatter(x=(i+1):(i+3), y=4:6 .* (10^i), yaxis="y$i", name="yaxis$i data")
    push!(traces, trace)
end

plot(
    traces,
    Layout(
        xaxis_domain=[0.3, 0.7],
        yaxis=attr(title="yaxis title", titlefont_color="red"),
        yaxis2=attr(
            title="yaxis2 title", titlefont_color="blue",
            overlaying="y", side="left", position=0.15, anchor="free"
        ),
        yaxis3=attr(
            title="yaxis3 title", titlefont_color="green",
            overlaying="y", side="right", anchor="x",
        ),
        yaxis4=attr(
            title="yaxis4 title", titlefont_color="orange",
            overlaying="y", side="right", position=0.85, anchor="free",
        ),
        title_text="multiple y-axes example",
    )
)