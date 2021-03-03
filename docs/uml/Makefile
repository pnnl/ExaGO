
all: opflow opflow-simple

opflow:
	cat opflow.mermaid opflow_defs.mermaid ps_defs.mermaid \
		| mmdc -o opflow-uml.png -H 3000 -w 3800

opflow-simple:
	cat opflow-simple.mermaid | mmdc -o opflow-simple-uml.png

clean:
	rm *.png
