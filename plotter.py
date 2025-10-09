
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("wavefunctions.csv")

plt.figure(figsize=(8,6))
for col in df.columns[1:]:
    plt.plot(df["x"], df[col], label=col)

plt.title("Quantum Harmonic Oscillator Wavefunctions")
plt.xlabel("x")
plt.ylabel("ψ(x)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
