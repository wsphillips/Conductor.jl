using Conductor, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, ms, mM, µM
import Conductor: Na, K, Ca, Cation, Leak

Vₘ = ParentScope(MembranePotential(-4.6mV))
Caᵢ = Concentration(Calcium, 0.2µM, dynamic = true)
ICa = IonCurrent(Calcium, aggregate = true)
#@parameters ϕ β = 0.075
#@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ]);

# Sodium channels
nav_kinetics = [
                # Activation
                Gate(AlphaBeta,
                     (0.32(13.1 - Vₘ)) / (exp((13.1 - Vₘ) / 4.0) - 1.0),
                     (0.28(Vₘ - 40.1)) / (exp((Vₘ - 40.1) / 5.0) - 1.0),
                     p = 2,
                     name = :m)
                # Inactivation
                Gate(AlphaBeta,
                     0.128 * exp((17.0 - Vₘ) / 18.0),
                     4.0 / (1.0 + exp((40.0 - Vₘ) / 5.0)),
                     name = :h)]

# Calcium channels
αᵣ = ifelse(Vₘ <= zero(Float64), 0.005, exp(-Vₘ / 20.0) / 200.0)

ca_kinetics = [
               # Activation
               Gate(AlphaBeta,
                    1.6 / (1.0 + exp(-0.072 * (Vₘ - 65.0))),
                    (0.02 * (Vₘ - 51.1)) / (exp((Vₘ - 51.1) / 5.0) - 1.0),
                    p = 2,
                    name = :s)
               # Inactivation; unused in Pinsk & Rinzel reduction
               Gate(AlphaBeta,
                    αᵣ,
                    ifelse(Vₘ <= zero(Float64), zero(Float64), 0.005 - αᵣ),
                    name = :r)]
# Delayed rectifier potassium
kdr_kinetics = [
    Gate(AlphaBeta,
         (0.016 * (35.1 - Vₘ)) / (exp((35.1 - Vₘ) / 5.0) - 1.0),
         0.25 * exp((20.0 - Vₘ) / 40.0),
         name = :n)]

# After-hyperpolarization potassium
# Note: The Traub model calls Calcium concentration 'χ' 
kahp_kinetics = [
    Gate(AlphaBeta,
         min(0.00002 * ParentScope(Caᵢ), 0.01),
         0.001,
         name = :q)]

αc = ifelse(Vₘ <= 50,
            (exp((Vₘ - 10.0) / 11.0) - exp((Vₘ - 6.5) / 27.0)) / 18.975,
            2.0 * exp((6.5 - Vₘ) / 27.0))

# Calcium-dependent potassium
kca_kinetics = [
    # Activation
    Gate(AlphaBeta,
         αc,
         ifelse(Vₘ <= 50,
                2.0 * exp((6.5 - Vₘ) / 27.0) - αc,
                zero(Float64)),
         name = :c),
    # Calcium saturation term
    Gate(SimpleGate, min(ParentScope(Caᵢ) / 250.0, 1.0),
         name = :χ)]

# A-type Potassium
ka_kinetics = [
               # Activation
               Gate(AlphaBeta,
                    (0.02 * (13.1 - Vₘ)) / (exp((13.1 - Vₘ) / 10.0) - 1.0),
                    (0.0175(Vₘ / 40.1)) / (exp((Vₘ - 40.1) / 10.0) - 1.0),
                    name = :a)
               # Inactivation
               Gate(AlphaBeta,
                    0.0016 * exp((-13.0 - Vₘ) / 18.0),
                    0.05 / (1 + exp((10.1 - Vₘ) / 5.0)),
                    name = :b)]

@named NaV = IonChannel(Sodium, nav_kinetics)
@named CaS = IonChannel(Calcium, ca_kinetics)
@named Kdr = IonChannel(Potassium, kdr_kinetics)
@named KAHP = IonChannel(Potassium, kahp_kinetics)
@named KA = IonChannel(Potassium, ka_kinetics)
@named KCa = IonChannel(Potassium, kca_kinetics)
@named leak = IonChannel(Leak)
