const { readFileSync, writeFileSync } = require('fs');
const { exec } = require('child_process');
const { fileSync } = require('tmp');

// function rscript(rfile, input) {
function rscript(rfile, expressionFile, genotypeFile, associationFile, gwasFile) {
    console.log('Files reached R-Wrapper:', expressionFile, genotypeFile, associationFile, gwasFile);
    return new Promise((resolve, reject) => {
        const workingDirectory = JSON.stringify(__dirname);
        console.log('R Working directory:', workingDirectory);

        expressionFile = JSON.stringify(expressionFile);
        genotypeFile = JSON.stringify(genotypeFile);
        associationFile = JSON.stringify(associationFile);
        gwasFile = JSON.stringify(gwasFile);
    
        var code = readFileSync(rfile).toString();
        
        // make sure the R statement below is not appended to a comment in R code file
        code += `eqtl(${workingDirectory}, ${genotypeFile}, ${expressionFile}, ${associationFile}, ${gwasFile})`;
        // console.log(code);

        const rcode = `
            suppressWarnings(suppressMessages(suppressPackageStartupMessages(
                jsonlite::toJSON({${code}}, auto_unbox=T)
            )))
        `;

        const tmpFile = fileSync();
        writeFileSync(tmpFile.name, rcode);

        const process = exec(
            `Rscript --vanilla "${tmpFile.name}"`, 
            { maxBuffer: 100 * 1024 * 1024 }, // increase default maxBuffer from ~4kb to 100mb
            (error, stdout, stderr) => {
                try {
                    if (stdout) {
                        resolve(JSON.parse(stdout.toString()));
                    } else {
                        if (error) reject(error);
                        if (stderr) reject(stderr);
                    }
                } catch(error) {
                    reject(error.toString());
                }
            }
        );
    });
}

module.exports = rscript;