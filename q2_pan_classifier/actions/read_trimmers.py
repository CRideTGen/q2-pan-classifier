def fastp_paired():
    cmd = f" fastp --in1 {inputs[0]} --in2 {inputs[1]}  --html {outputs[2]} --json {outputs[3]} --out1 {outputs[0]} --out2 {outputs[1]}; "
