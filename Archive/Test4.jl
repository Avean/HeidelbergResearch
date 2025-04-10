function mandel(c)
    z = c
    maxiter = 80
    for n in 1:maxiter
    if abs(z) > 2
    return n - 1
    end
    z = z^2 + c
    end
    return maxiter
   end
   
   function randmsq()
    x = rand(10000, 1000)
    y = mean(x.^2, dims=1)
    return y
   end

#    function randmsq_timed()
#     @timeit to "randmsq" begin
#     x = @timeit to "rand" rand(10000, 1000)
#     y = @timeit to "mean" mean(x.^2, dims=1)
#     return y
#     end
#     end


function create_matrix_struct(original_struct::Type{T}) where T
    # Pobieramy nazwy pól oryginalnej struktury
    fields = fieldnames(T)
    
    # Pobieramy typy pól oryginalnej struktury
    types = fieldtypes(T)
    
    # Generujemy nazwę nowej struktury
    newname = Symbol(nameof(T), :Matrix)
    
    # Przygotowujemy wyrażenie dla nowej struktury
    struct_expr = Expr(:struct, false, newname, 
        Expr(:block, 
            [Expr(:(::), field, Expr(:curly, :Matrix, type)) for (field, type) in zip(fields, types)]...
        )
    )
    print(struct_expr)
    # Ewaluujemy dynamicznie wygenerowaną strukturę
    eval(struct_expr)
    
    # Zwracamy typ nowo utworzonej struktury
    return eval(newname)
end

# Przykładowa oryginalna struktura
struct OriginalStruct
    field1::Int
    field2::Float64
    field3::String
end

# Utworzenie struktury z macierzami
MatrixStruct = create_matrix_struct(OriginalStruct)

# Demonstracja użycia
matrix_instance = MatrixStruct(
    [1 2; 3 4],           # field1 jako Matrix{Int}
    [1.0 2.0; 3.0 4.0],   # field2 jako Matrix{Float64}
    ["a" "b"; "c" "d"]    # field3 jako Matrix{String}
)

# Wyświetlenie struktury
println(matrix_instance)

# Dostęp do konkretnych pól
println(matrix_instance.field1)
println(matrix_instance.field2)
println(matrix_instance.field3)

# a = 2;
# b = 1;
# c = 4;
# z1 = :($a+$b+$c)

# z = quote a+b+c end

# typeof(z.args)

# macro twostep(arg)
#     println("I execute at parse time. The argument is: ", arg)
#     return :(println("I execute at runtime. The argument is: ", $arg))
# end

# @twostep (macro with 1 method)



            # input_struct = params
            # # Pobieramy oryginalną strukturę i jej nazwy pól
            # original_type = typeof(input_struct)
            # field_names = fieldnames(original_type)
            
            # # Dynamiczne tworzenie nazwy nowej struktury
            # new_struct_name = Symbol(string(original_type) * "Matrix")
            
            # # Dynamiczne tworzenie definicji nowej struktury z polami Matrix{Float64}
            # struct_expr = Expr(:struct, false, new_struct_name, 
            #     Expr(:block, 
            #         map(field_names) do field_name
            #             Expr(:(::), field_name, :(Matrix{Float64}))
            #         end...
            #     )
            # )
            
            # # Wykonanie definicji struktury
            # W = eval(struct_expr)
            
            # # Przygotowanie wartości dla nowej struktury
            # field_values = map(field_names) do field_name
            #     original_value = getfield(input_struct, field_name)
            #     return Matrix{Float64}(fill(Float64(original_value), 1, 1))
            # end
            
            # # Dynamiczne utworzenie instancji nowej struktury
            # new_struct_type = eval(new_struct_name)
            # return new_struct_type(field_values...)

function convert_struct_to_matrix(input_struct)
    # Pobieramy oryginalną strukturę i jej nazwy pól
    original_type = typeof(input_struct)
    field_names = fieldnames(original_type)
    
    # Dynamiczne tworzenie nazwy nowej struktury
    new_struct_name = Symbol(string(original_type) * "Matrix")
    
    # Dynamiczne tworzenie definicji nowej struktury z polami Matrix{Float64}
    struct_expr = Expr(:struct, false, new_struct_name, 
        Expr(:block, 
            map(field_names) do field_name
                Expr(:(::), field_name, :(Matrix{Float64}))
            end...
        )
    )
    
    # Wykonanie definicji struktury
    eval(struct_expr)
    
    # Przygotowanie wartości dla nowej struktury
    field_values = map(field_names) do field_name
        original_value = getfield(input_struct, field_name)
        return Matrix{Float64}(fill(Float64(original_value), 1, 1))
    end
    
    # Dynamiczne utworzenie instancji nowej struktury
    new_struct_type = eval(new_struct_name)
    return new_struct_type(field_values...)
end
            
    struct Parameters2
        D::Float64
        κ::Float64
        A::Int
    end
    
    params = Parameters2(0.1, 10.0, 15)
    matrix_params = convert_struct_to_matrix(params)

    matrix_params2 = convert_struct_to_matrix(matrix_params)
   