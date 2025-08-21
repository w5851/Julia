"""
PNJL 物理仿真软件包

一个用于研究QCD相变的综合Julia软件包。
"""
module PNJLPhysicsSimulation

# 核心模块
include("core/constants.jl")
include("core/math_utils.jl")
include("core/integration_interface.jl")
include("core/thermodynamics.jl")
# include("core/function_registry.jl")  # 函数注册表（已移除）
include("core/autodiff_interface.jl")  # 专门的自动微分模块（已重命名）

# 模型模块 - 必须在统一接口之前包含
include("models/gas_liquid/constants.jl")
include("models/gas_liquid/functions.jl")
include("models/pnjl/constants.jl")
include("models/pnjl/functions.jl")
include("models/pnjl_aniso/constants.jl")
include("models/pnjl_aniso/functions.jl")
include("models/rotation/constants.jl")
include("models/rotation/functions.jl")

# 高层接口 - 在模型之后包含（部分接口已精简，相关实现移至独立模块）
# include("core/model_configuration.jl")
# include("core/unified_physics_interface.jl")  # 统一物理接口
# 以下公共接口文件目前仅作为存储文件保留，暂不在包加载时包含，
# 以避免在预编译阶段触发未完成的接口或依赖引发错误。
# 如需启用，请确认接口稳定并取消下面的注释：
# include("core/unified_physics_public_interface.jl")  # 新的通用公共接口（已注释）

# 引入并 re-export 独立的求解器模块（提供 solve_equilibrium_equations）
include("core/solver_interface.jl")

# 导出核心功能和接口
export PhysicalConstants, IntegrationInterface, ModelConfiguration, 
       Thermodynamics, MathUtils,
       GasLiquidConstants, GasLiquidFunctions,
       PNJLConstants, PNJLFunctions, 
       PNJLAnisoConstants, PNJLAnisoFunctions,
       RotationConstants, RotationFunctions,
       # UnifiedPhysicsInterface, (已移除)
       solve_equilibrium_equations

end  # module PNJLPhysicsSimulation
