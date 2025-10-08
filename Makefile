build:
    latexmk -pdf -output-directory=build main
clean:
    latexmk -C -output-directory=build