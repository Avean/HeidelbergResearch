
input_struct = Sets.Set1.Diff

LapMat = Lap.Dir2
EigVec = Eig.S
dt = 0.1;


    field_names = fieldnames(typeof(input_struct))
    
    new_struct_name = Symbol("DiffMat2")
    

    struct_expr = Expr(:struct, false, new_struct_name, 
        Expr(:block, 
            map(field_names) do field_name
                Expr(:(::), field_name, :DiffusionMat)   
            end...
        )
    )
    display(struct_expr)

    eval(struct_expr)
            
    # Nowy typ struktury
   new_struct_type = eval(new_struct_name)
   return (map(h -> FillDiffMatrix(getfield(input_struct,h), LapMat, EigVec, dt),field_names));


