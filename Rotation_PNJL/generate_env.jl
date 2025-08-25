using Pkg, UUIDs, TOML

PROJECT_DIR = @__DIR__

# 1. 初始化项目
Pkg.activate(PROJECT_DIR)

# 创建 Project.toml 文件
project_data = Dict(
    "name" => basename(PROJECT_DIR),
    "uuid" => string(uuid1()),
    "version" => "0.1.0"
)

using TOML
TOML.tryparsefile(joinpath(PROJECT_DIR, "Project.toml")) == nothing && 
    TOML.writefile(joinpath(PROJECT_DIR, "Project.toml"), project_data)

# 2. 分析 src 目录中的 Julia 文件并添加常见依赖
if isdir(joinpath(PROJECT_DIR, "src"))
    println("Found src directory, analyzing Julia files...")
    
    # 读取所有 .jl 文件并查找 using 和 import 语句
    julia_files = filter(f -> endswith(f, ".jl"), readdir(joinpath(PROJECT_DIR, "src"), join=true))
    
    dependencies = Set{String}()
    for file in julia_files
        try
            content = read(file, String)
            # 查找 using 和 import 语句
            for line in split(content, '\n')
                line = strip(line)
                if startswith(line, "using ") || startswith(line, "import ")
                    # 提取包名
                    parts = split(replace(line, r"^(using|import)\s+" => ""), r"[,\s:]")
                    for part in parts
                        pkg = strip(part)
                        if !isempty(pkg) && !startswith(pkg, ".")
                            push!(dependencies, pkg)
                        end
                    end
                end
            end
        catch e
            println("Warning: Could not read file $file: $e")
        end
    end
    
    # 添加发现的依赖
    for dep in dependencies
        if !in(dep, ["Base", "Core", "Main"])
            try
                println("Adding dependency: $dep")
                Pkg.add(dep)
            catch e
                println("Warning: Could not add dependency $dep: $e")
            end
        end
    end
end

# 3. 添加常见元数据
open(joinpath(PROJECT_DIR, "Project.toml"), "a") do io
    write(io, """
    authors = ["Your Name <you@example.com>"]
    
    [compat]
    julia = "^$(VERSION.major).$(VERSION.minor)"  # 自动获取当前Julia版本
    """)
end

println("Environment created at $PROJECT_DIR")