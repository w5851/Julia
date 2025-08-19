"""
PNJL 物理仿真软件包

一个用于研究QCD相变的综合Julia软件包。
"""
module PNJLPhysicsSimulation

# 核心模块
include("core/constants.jl")
include("core/math_utils.jl")
include("core/integration.jl")
include("core/integration_interface.jl")
include("core/thermodynamics.jl")
include("core/function_registry.jl")  # 函数注册表

# 模型模块 - 必须在统一接口之前包含
include("models/gas_liquid/constants.jl")
include("models/gas_liquid/functions.jl")
include("models/pnjl/constants.jl")
include("models/pnjl/functions.jl")
include("models/pnjl_aniso/constants.jl")
include("models/pnjl_aniso/functions.jl")
include("models/rotation/constants.jl")
include("models/rotation/functions.jl")

# 高层接口 - 在模型之后包含
include("core/model_configuration.jl")
# include("core/unified_physics_interface.jl")  # 暂时注释掉避免依赖问题

# 暂时导出所有内容进行测试
export PhysicalConstants, Integration, IntegrationInterface, ModelConfiguration, 
       FunctionRegistry, Thermodynamics, MathUtils,
       GasLiquidConstants, GasLiquidFunctions,
       PNJLConstants, PNJLFunctions, 
       PNJLAnisoConstants, PNJLAnisoFunctions,
       RotationConstants, RotationFunctions

end  # module PNJLPhysicsSimulation
