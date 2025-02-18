using LinearAlgebra
using DataFrames
using CSV
using Plots
using SparseArrays

# Calcular la matriz de admitancia nodal

function calcular_ybus(lines,nodes)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    Ybus : Matriz
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Ybus = zeros(num_nodes, num_nodes)*1im

    for k = 1:num_lines
        # Nodo de envío
        n1 = lines.FROM[k]
        # Nodo de recibo
        n2 = lines.TO[k]
        # Admitancia de la línea
        yL = 1/(lines.R[k]+lines.X[k]*1im)
        # Susceptancia de la línea
        Bs = lines.B[k]*1im/2
        # Valor del TAP
        t = lines.TAP[k]
        if lines.TAP[k] == 0
            Ybus[n1,n1] += yL + Bs   # Dentro de la diagonal
            Ybus[n1,n2] -= yL        # Fuera de la diagonal
            Ybus[n2,n1] -= yL        # Fuera de la diagonal
            Ybus[n2,n2] += yL + Bs   # Dentro de la diagonal
        else
            Ybus[n1,n1] += (t^2 - t)*yL  # Dentro de la diagonal
            Ybus[n1,n2] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n1] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n2] += (1-t)*yL      # Dentro de la diagonal
        end
    end
    return Ybus
end

# Se calcula la matriz Bbus
function calcular_bbus(lines,nodes)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
    Salida :    Bbus : Matriz
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Bbus = zeros(num_nodes, num_nodes)
    for k = 1:num_lines
        # Nodo de envío
        n1 = lines.FROM[k]
        # Nodo de recibo
        n2 = lines.TO[k]
        # Admitancia de la línea
        BL = 1/(lines.X[k])
        Bbus[n1,n1] += BL        # Dentro de la diagonal
        Bbus[n1,n2] -= BL        # Fuera de la diagonal
        Bbus[n2,n1] -= BL        # Fuera de la diagonal
        Bbus[n2,n2] += BL        # Dentro de la diagonal
    end
    return Bbus
end

function flujo_dc(lines,nodes, Bbus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
                Bbus  : Matriz 
    Salida :    df_theta : DataFrame de los angulos nodales
                df_pkm : DataFrame del flujo de potencia por las líneas
    """
    s = nodes[nodes.TYPE .== 3, "NUMBER"]
    Bbus = Bbus[setdiff(1:end, s), setdiff(1:end, s)]
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    # Se calculan las Potencias nodales.
    Pn = zeros(num_nodes)
    for k = 1:num_nodes
        Pn[k]= nodes.PGEN[k]-nodes.PLOAD[k]
    end
    Pn = Pn[setdiff(1:end, s)]

    # Calculo de los angulos nodales.
    thetas_sin_slack = inv(Bbus)*Pn

    # Se implementa el nodo Slack
    theta = zeros(nrow(nodes))
    theta[1:end .!= s] = thetas_sin_slack

    # Calcular los Pkm.
    Pkm = zeros(num_lines)
    for k = 1:num_lines
        Pkm[k] = (theta[lines.FROM[k]]-theta[lines.TO[k]])/lines.X[k]
    end

    df_theta = DataFrame(Nodo = nodes.NUMBER, Theta = theta)
    df_pkm = DataFrame(From = lines.FROM, To = lines.TO, Pkm = Pkm)

    return df_theta, df_pkm
end

function flujo_dc_data(lines,nodes, Bbus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
                Bbus  : Matriz 
    Salida :    df_theta : DataFrame de los angulos nodales
                df_pkm : DataFrame del flujo de potencia por las líneas
    """
    s = nodes[nodes.TYPE .== 3, "NUMBER"]
    Bbus = Bbus[setdiff(1:end, s), setdiff(1:end, s)]
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    # Se calculan las Potencias nodales.
    Pn = zeros(num_nodes)
    for k = 1:num_nodes
        Pn[k]= nodes.PGEN[k]-nodes.PLOAD[k]
    end
    Pn = Pn[setdiff(1:end, s)]

    # Calculo de los angulos nodales.
    thetas_sin_slack = inv(Bbus)*Pn

    # Se implementa el nodo Slack
    theta = zeros(nrow(nodes))
    theta[1:end .!= s] = thetas_sin_slack

    # Calcular los Pkm.
    Pkm = zeros(num_lines)
    for k = 1:num_lines
        Pkm[k] = (theta[lines.FROM[k]]-theta[lines.TO[k]])/lines.X[k]
    end
    return theta, Pkm
end


# Función principal
lines = DataFrame(CSV.File("lines.csv"))
nodes = DataFrame(CSV.File("nodes.csv"))

### Cálculo de la Matriz Ybus
Ybus = calcular_ybus(lines,nodes)
Bbus = calcular_bbus(lines,nodes)

### Cálculo del flujo DC
(theta,pkm) = flujo_dc(lines,nodes,Bbus) 
display(theta)
display(pkm)

### Analisis de contingencias N-1
function contingencias(lines,nodes, Bbus)
    """
    Entradas:   lines: DataFrames
                nodes : DataFrames
                Bbus  : Matriz 
    Salida :    theta_cont : DataFrame de los angulos nodales para cada caso evaluado.
                pkm_cont : DataFrame del flujo de potencia por las líneas para cada caso evaluado.
    """
    num_lines = nrow(lines)
    num_nodes = nrow(nodes)
    theta_cont = DataFrame()  
    pkm_cont = DataFrame()
    for k = 1:num_lines
        Bbus_temp = copy(Bbus)
        # Se modifica la matriz Bbus sacando una de las líneas
        Bbus_temp[lines.FROM[k],lines.TO[k]] += 1/lines.X[k] # fuera de la diagonal
        Bbus_temp[lines.TO[k],lines.FROM[k]] += 1/lines.X[k] # fuera de la diagonal
        Bbus_temp[lines.FROM[k],lines.FROM[k]] -= 1/lines.X[k] # dentro de la diagonal
        Bbus_temp[lines.TO[k],lines.TO[k]] -= 1/lines.X[k] # dentro de la diagonal
        # Se ejecuta el flujo de carga
        Bbus_det = copy(Bbus_temp)
        s = nodes[nodes.TYPE .== 3, "NUMBER"]
        Bbus_det = Bbus_det[setdiff(1:end, s), setdiff(1:end, s)]
        if det(Bbus_det) == 0
            theta = zeros(num_nodes)
            pkm = zeros(num_lines)
            theta[1:end] .= Inf64
            pkm[1:end] .= Inf64
        else
            theta,pkm = flujo_dc_data(lines,nodes,Bbus_temp) 
        end
        # Se modifica la posición k debido a que en esta la linea no tiene flujo de potencia
        pkm[k] = 0

        # Se condiciona el guardado de cada uno de los casos debido a que el sistema puede poseer lineas en paralelo
        if k != 1 && (lines.FROM[k] == lines.FROM[k-1]) && (lines.TO[k] == lines.TO[k-1])
            col_theta = "Angulos caso línea $(lines.FROM[k]) - $(lines.TO[k]) paralela"
            theta_cont[!,col_theta] = theta
            col_pkm = "Flujos caso línea $(lines.FROM[k]) - $(lines.TO[k]) paralela"
            pkm_cont[!,col_pkm] = pkm
        else
            col_theta = "Angulos caso línea $(lines.FROM[k]) - $(lines.TO[k])"
            theta_cont[!,col_theta] = theta
            col_pkm = "Flujos caso línea $(lines.FROM[k]) - $(lines.TO[k])"
            pkm_cont[!,col_pkm] = pkm
        end
    end
    return theta_cont, pkm_cont
end

theta_cont, pkm_cont = contingencias(lines,nodes,Bbus)
CSV.write("angulos_nodales.csv", theta_cont)
CSV.write("flujo_de_potencia.csv", pkm_cont)


### Analisis gráfico de los resultados obtenidos.

df_datatheta = CSV.read("angulos_nodales.csv", DataFrame)
df_datapkm =  CSV.read("flujo_de_potencia.csv", DataFrame)

# Convertir el DataFrame a una matriz
matriz_theta = Matrix(df_datatheta)
matriz_pkm = Matrix(df_datapkm)


display(matriz_pkm)

# Creando el mapa de calor para los angulos
p1 = heatmap(
    matriz_theta,
    c=:plasma,
    aspect_ratio=:equal,
    xlabel="Casos",
    ylabel="Nodo",
    title="Angulos para las contingencias n-1"
)

# Creando el mapa de calor para los flujos
p2 = heatmap(
    matriz_pkm,
    c=:plasma,
    aspect_ratio=:equal,
    xlabel="Casos",
    ylabel="Lineas",
    title="Flujos de potencia para las lineas"
)

# Combinar los plots
plot(p1, p2, layout=(1,2), size=(1200,500))
savefig("Mapa_de_calor.png")
