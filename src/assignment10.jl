#!/usr/bin/env julia
#
# Copyright 2022 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
module assignment10

function z(P, params)
    #Kareem, L.A., Iwalewa, T.M. & Al-Marhoun, M. New explicit correlation for the 
    #compressibility factor of natural gas: linearized z-factor isotherms. J Petrol 
    #Explor Prod Technol 6, 481–492 (2016). https://doi.org/10.1007/s13202-015-0209-3

    ρ, T = params
    
    γg = ρ / 0.0806
    
    Tpc = (169.2 + 349.5 * γg - 74.0 * γg ^ 2) 
    Ppc = (756.8 - 131.07* γg - 3.6 * γg ^ 2) 
    
    Tpr = T / Tpc
    Ppr = P / Ppc
    
    t = 1 / Tpr
    
    a₁  = 0.317842
    a₂  = 0.382216
    a₃  = -7.768354
    a₄  = 14.290531
    a₅  = 0.000002
    a₆  = -0.004693
    a₇  = 0.096254
    a₈  = 0.166720
    a₉  = 0.966910
    a₁₀ = 0.063069
    a₁₁ = -1.966847
    a₁₂ = 21.0581
    a₁₃ = -27.0246
    a₁₄ = 16.23
    a₁₅ = 207.783
    a₁₆ = -488.161
    a₁₇ = 176.29
    a₁₈ = 1.88453
    a₁₉ = 3.05921
    
    A = a₁ * t * exp(a₂ * (1 - t) ^ 2) * Ppr
    B = a₃ * t + a₄ * t ^ 2 + a₅ * t ^ 6 * Ppr ^ 6
    C = a₉ + a₈ * t * Ppr + a₇ * t ^ 2 * Ppr ^ 2 + a₆ * t ^ 3 * Ppr ^ 3
    D = a₁₀ * t * exp(a₁₁ * (1 - t) ^ 2)
    E = a₁₂ * t + a₁₃ * t ^ 2 + a₁₄ * t ^ 3
    F = a₁₅ * t + a₁₆ * t ^ 2 + a₁₇ * t ^ 3
    G = a₁₈ + a₁₉ * t
    
    y = D * Ppr / ((1 + A ^ 2) / C - A ^ 2 * B / C ^ 3)
    
    D * Ppr * (1 + y + y ^ 2 - y ^ 3) / (D * Ppr + E * y ^ 2 - F * y ^ G) / (1 - y) ^ 3
end


function gas_solver(Pₒ, Vₚ, Pₛ, Tₛ, Pwf,  n, J, T, ρ, tmax)
    #return the solver composite type (what solve() returns)
end

export gas_solver

end
