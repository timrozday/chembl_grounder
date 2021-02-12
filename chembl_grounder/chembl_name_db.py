import sqlalchemy as sa
from sqlalchemy.ext import declarative
import sqlalchemy.ext.associationproxy as sa_ap
from sqlalchemy.orm import foreign, remote

class ChemblNameDB():
    Base = sa.ext.declarative.declarative_base()
    
    def __init__(self, db_path=""):
        self.db_path = db_path
        self.engine = sa.create_engine(self.db_path)  # pool_recycle=3600
        self.SessionMaker = sa.orm.sessionmaker(bind=self.engine)
        
    def create(self):
        self.Base.metadata.create_all(self.engine)
        
    def drop_all(self):
        self.Base.metadata.drop_all(self.engine)
    
    class Substance(Base):
        __tablename__ = "substance"
        
        id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
        chembl_id = sa.Column(sa.String(10), index=True)
        molregno = sa.Column(sa.Integer, index=True)
        drugbase_id = sa.Column(sa.Integer, index=True)
        
        name_assocs = sa.orm.relationship("Substance2Name", back_populates="substance", cascade="save-update, merge, delete")
        
        @classmethod
        def unique(cls, session, chembl_id=None, molregno=None, drugbase_id=None, flush=True):
            obj = cls.get_unique(session, chembl_id=chembl_id, molregno=molregno, drugbase_id=drugbase_id, flush=False)
            if chembl_id:
                obj.chembl_id = chembl_id
            if molregno:
                obj.molregno = molregno
            if drugbase_id:
                obj.drugbase_id = drugbase_id
            
            if flush:
                session.flush()

            return obj
        
        @classmethod
        def get_unique(cls, session, chembl_id=None, molregno=None, drugbase_id=None, flush=True):
            if drugbase_id:
                obj = session.query(cls).filter(cls.drugbase_id==drugbase_id).first()
                if obj:
                    return obj

            if chembl_id:
                obj = session.query(cls).filter(cls.chembl_id==chembl_id).first()
                if obj:
                    return obj

            if molregno:
                obj = session.query(cls).filter(cls.molregno==molregno).first()
                if obj:
                    return obj

            obj = cls(chembl_id=chembl_id, molregno=molregno, drugbase_id=drugbase_id)
            session.add(obj)
            if flush:
                session.flush()

            return obj
    
    class Substance2Name(Base):
        __tablename__ = "substance2name"
        __table_args__ = (sa.Index('substance2name_substanceid_nameid_idx', 'substance_id', "name_id", unique=True), )
        
        substance_id = sa.Column(sa.Integer, sa.ForeignKey("substance.id"), primary_key=True)
        name_id = sa.Column(sa.Integer, sa.ForeignKey("name.id"), primary_key=True)
        
        substance = sa.orm.relationship("Substance", back_populates="name_assocs", cascade="save-update, merge")
        name = sa.orm.relationship("Name", back_populates="substance_assocs", cascade="save-update, merge")
        
        @classmethod
        def unique(cls,session,flush=True,**kw):
            return ChemblNameDB._unique(session, \
                                        cls, \
                                        flush=flush, \
                                        hashfunc=lambda substance_id,name_id:f'Substance2Name:{substance_id};{name_id}', \
                                        queryfunc=lambda q,substance_id,name_id:q.filter(cls.substance_id==substance_id,\
                                                                                         cls.name_id==name_id), \
                                        kw=kw)
    
    class Name(Base):
        __tablename__ = "name"
        __table_args__ = (sa.Index('name_name_type_idx', 'name', 'table', "type", unique=True), )
        
        id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
        name = sa.Column(sa.String(100), index=True)
        lower = sa.Column(sa.String(100), index=True)
        table = sa.Column(sa.String(16), index=True)
        type = sa.Column(sa.String(16), index=True)
        
        substance_assocs = sa.orm.relationship("Substance2Name", back_populates="name", cascade="save-update, merge, delete")
        
        @classmethod
        def unique(cls,session,flush=True,**kw):
            return ChemblNameDB._unique(session, \
                                        cls, \
                                        flush=flush, \
                                        hashfunc=lambda name,table,type:f'Name:{name};{table};{type}', \
                                        queryfunc=lambda q,name,table,type:q.filter(cls.name==name,\
                                                                              cls.table==table,\
                                                                              cls.type==type), \
                                        kw=kw)
    
    def _unique(session, cls, flush=True, hashfunc=lambda name:name, queryfunc=lambda query,name:query.filter(Table.name==name), arg=[], kw={}, constructor=None):
        if not constructor:
            constructor = cls

        cache = getattr(session, '_unique_cache', None)
        if cache is None:
            session._unique_cache = cache = {}

        key = (cls, hashfunc(*arg, **kw))
        if key in cache:
            return cache[key]
        else:
            with session.no_autoflush:
                q = session.query(cls)
                q = queryfunc(q, *arg, **kw)
                obj = q.first()
                if not obj:
                    obj = constructor(*arg, **kw)
                    session.add(obj)
                    if flush:
                        session.flush()
            
            cache[key] = obj
            return obj
