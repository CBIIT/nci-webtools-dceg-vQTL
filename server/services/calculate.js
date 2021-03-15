const { v1: uuidv1 } = require('uuid');
const r = require('r-wrapper').async;
const path = require('path');
const logger = require('../services/logger');

async function qtlsCalculateMain(params, res, next) {
    logger.info("Execute /qtls-calculate-main");
    logger.debug(params);
    const {
        select_qtls_samples,
        select_gwas_sample,
        associationFile,
        quantificationFile,
        genotypeFile,
        gwasFile,
        LDFile,
        select_pop,
        select_gene,
        select_dist,
        select_ref,
        recalculateAttempt,
        recalculatePop,
        recalculateGene,
        recalculateDist,
        recalculateRef,
        workingDirectory
    } = params;
    const rfile = path.resolve(__dirname, 'query_scripts', 'QTLs', 'qtls.r');
    const request = uuidv1();
    try {
        const wrapper = await r(
            path.resolve(__dirname, 'query_scripts', 'wrapper.R'),
            "qtlsCalculateMain",
            [
                rfile,
                workingDirectory,
                select_qtls_samples.toString(),
                select_gwas_sample.toString(),
                associationFile,
                quantificationFile,
                genotypeFile,
                gwasFile,
                LDFile,
                request,
                select_pop.toString(),
                select_gene,
                select_dist,
                select_ref.toString(),
                recalculateAttempt,
                recalculatePop,
                recalculateGene,
                recalculateDist,
                recalculateRef
            ]
        );
        logger.info("Finished /qtls-calculate-main");
        res.json(JSON.parse(wrapper));
    } catch (err) {
        logger.error("Error /qtls-calculate-main");
        logger.error(err);
        res.status(500).json(err);
    }
}

module.exports = {
    qtlsCalculateMain
}