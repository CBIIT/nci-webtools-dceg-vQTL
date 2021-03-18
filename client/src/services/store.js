import { createStore, applyMiddleware, compose } from 'redux';
import { rootReducer } from './reducers';
import ReduxThunk from 'redux-thunk';

export const getInitialState = async () => {
  let initialState = {
    qtlsGWAS: {
      openSidebar: true,
      select_qtls_samples: false,
      select_gwas_sample: false,
      associationFile: null,
      quantificationFile: null,
      genotypeFile: null,
      gwasFile: null,
      LDFile: null,
      request: '',
      select_pop: false,
      select_gene: false,
      select_dist: 100,
      select_ref: false,
      recalculateAttempt: false,
      recalculatePop: false,
      recalculateGene: false,
      recalculateDist: false,
      recalculateRef: false,
      submitted: false,
      isLoading: false,
      isError: false,
      activeResultsTab: 'locus-alignment',
      top_gene_variants: null,
      all_gene_variants: null,
      gene_list: null,
      inputs: null,
      messages: null,
      locus_quantification: null,
      locus_quantification_heatmap: null,
      locus_alignment: {
        data: null,
        layout: null,
        top : {
          chr: null,
          gene_symbol: null
        }
      },
      locus_alignment_gwas_scatter: null,
      locus_colocalization_correlation: null,
      gwas: null,
      locus_table: null
      // locus_table: {
      //   globalFilter: '',
      //   data: [],
      //   pagination: {
      //     pageIndex: 0,
      //     pageSize: 0
      //   }
      // }
    },
    errorModal: {
      visible: false,
      details: ``,
      message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:FORGE2-TFWebAdmin@cancer.gov">ezQTLWebAdmin@cancer.gov</a>.`,
    },
  };

  return initialState;
}

export const getStore = async _ => createStore(
  rootReducer,
  await getInitialState(),
  compose(
    applyMiddleware(ReduxThunk),
    window.__REDUX_DEVTOOLS_EXTENSION__
      ? window.__REDUX_DEVTOOLS_EXTENSION__({ trace: true })
      : e => e
  )
);