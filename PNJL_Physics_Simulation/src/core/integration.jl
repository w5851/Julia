"""
物理仿真核心积分工具。

此模块为PNJL仿真包中的不同物理模型提供
通用积分函数。
"""
module Integration

using FastGaussQuadrature

export gauleg

"""
    gauleg(a, b, n) -> (nodes, weights)

高斯-勒让德求积积分节点和权重。

# 参数
- `a::Real`: 积分下界
- `b::Real`: 积分上界  
- `n::Int`: 积分点数

# 返回值
- `nodes::Vector{Float64}`: 积分节点
- `weights::Vector{Float64}`: 积分权重

# 示例
```julia
nodes, weights = gauleg(0.0, 1.0, 10)
```
"""
function gauleg(a, b, n)
    t_nodes, t_weights = gausslegendre(n)
    nodes   = @. (b - a)/2 * t_nodes + (a + b)/2
    weights = @. (b - a)/2 * t_weights
    return nodes, weights
end

end  # module Integration
