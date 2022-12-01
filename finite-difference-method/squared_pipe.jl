using Plots

hb = 28.
ha = 350.
k = 2.1
Dx = 3.
Dy = 3.
Te = 30.
Ti = 220.

T = zeros(29, 29)

function squared_pipe!(T)

    #pontas externas
    T[1, 29] = (T[1, 28] + T[2, 29] + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+1))
    T[1, 1] = (T[2, 1] + T[1, 2] + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+1))
    T[29, 1] = (T[28, 1] + T[29, 2] + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+1))
    T[29, 29] = (T[29, 28] + T[28, 29] + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+1))

    #pontas internas
    T[6, 6] = (T[7, 6] + T[6, 7] + 2*(T[5, 6] + T[6, 5]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+3))
    T[6, 24] = (T[7, 24] + T[6, 23] + 2*(T[5, 24] + T[6, 25]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+3))
    T[24, 24] = (T[23, 24] + T[24, 23] + 2*(T[5, 6] + T[6, 5]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+3))
    T[24, 6] = (T[23, 6] + T[24, 7] + 2*(T[24, 5] + T[25, 6]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+3))

    #paredes com convecção

    #vertical esquerda externa
    for i in 1:27
        T[i+1, 1] = (T[i, 1] + T[i+2, 1] + 2*(T[i+1, 2]) + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+2))
    end

    #vertical externa direita
    for j in 1:27
        T[j+1, 29] = (T[j, 29] + T[j+2, 29] + 2*(T[j+1, 28]) + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+2))
    end

    #parede interna esquerda
    for i in 1:17
        T[i+6, 6] = (T[i+5, 6] + T[i+7, 6] + 2*(T[i+6, 5]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+2))
    end

    #parede interna direita
    for j in 1:17
        T[j+6, 24] = (T[j+5, 24] + T[j+7, 24] + 2*(T[j+6, 25]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+2))
    end

    #horizontal cima extern0
    for i in 1:27
        T[1, i+1] = (T[1, i] + T[1, i+2] + 2*(T[2, i+1]) + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+2))
    end

    #horizontal baixo externo
    for j in 1:27
        T[29, j+1] = (T[29, j] + T[29, j+2] + 2*(T[28, j+1]) + (2*hb*Dx*Te)/k)/(2*(((hb*Dx)/k)+2))
    end

    #horizontal cima interna
    for i in 1:17
        T[6, i+6] = (T[6, i+5] + T[6, i+7] + 2*(T[5, i+6]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+2))
    end

    #horizontal baixo interna
    for i in 1:17
        T[24, i+6] = (T[24, i+5] + T[24, i+7] + 2*(T[25, i+6]) + (2*ha*Dx*Ti)/k)/(2*(((ha*Dx)/k)+2))
    end

    #bloco interno cima
    for i in 1:4
        for j in 1:27
            T[i+1, j+1] = (T[i, j+1] + T[i+2, j+1] + T[i+1, j] + T[i+1, j+2])/4
        end
    end

    #bloco interno de baixo
    for i in 1:4
        for j in 1:27
            T[i+24, j+1] = (T[i+23, j+1] + T[i+25, j+1] + T[i+24, j] + T[i+24, j+2])/4
        end
    end

    #bloco lateral esquerdo
    for i in 1:19
        for j in 1:4
            T[i+5, j+1] = (T[i+4, j+1] + T[i+6, j+1] + T[i+5, j] + T[i+5, j+2])/4
        end
    end

    #bloco lateral direito
    for i in 1:19
        for j in 1:4
            T[i+5, j+24] = (T[i+4, j+24] + T[i+6, j+24] + T[i+5, j+23] + T[i+5, j+25])/4
        end
    end

    #return T
end

for i in 1:100
    resultado = squared_pipe!(T)
end

S = reverse(T, dims=1)

texts = [(j, i, text(round(S[i, j], digits=2), 2, :black, :center))
         for i in 1:29 for j in 1:29]

heatmap(T, aspect_ratio=1.0, color=:thermal, xlim=(0, 30))
annotate!(texts, linecolor=:blue)
savefig("pipe_heatmap.pdf")